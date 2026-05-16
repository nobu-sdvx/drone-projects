"""
propeller_design_analysis.py
============================

§5 プロペラ設計 — メイン解析(空力性能)

対象: 5/18 第2回審査会の概略計算書 §5
作成: 2026-05-16
設計根拠: Propeller_design_spec_v0.3 / user-request_propeller.md §3.1

役割:
  自作プロペラ(D=36mm, P=27mm, Crazyflie 45mm 形状スケーリング)の
  空力性能を 2 つの独立手段で評価し交差検証する:
    (1) QPROP 1.22(MIT Drela、BEM ソルバ)   ← 主たる数値
    (2) 自作 BEM(ブレード要素運動量理論)    ← 念のための検算(±20% 想定)
  併せて XFLR5 ポーラでブレード断面翼型 4 候補を動作 Re 域で比較し 1 つに確定する。

入力:
  xfoil_data/propeller_blade/*_T1_Re*.txt   … ブレード断面翼型ポーラ(XFLR5)
  qprop_data/powerup_propeller.prop          … QPROP プロペラ定義(暫定翼型)
  qprop_data/powerup_motor.mot               … QPROP モーター定義
  uiuc_propeller_data/.../cfnq_45p1_geom.txt … Crazyflie 形状(出典)

出力(コンソール + plots/):
  - propeller_3view.png            プロペラ三面図
  - blade_chord_distribution.png   弦長分布 c(r)(自作 vs Crazyflie)
  - blade_twist_distribution.png   ねじり角分布 β(r)(自作 vs Crazyflie)
  - airfoil_polars.png             ブレード断面翼型ポーラ(XFLR5)
  - thrust_vs_rpm_curves.png       推力 vs RPM(V=0,5,7.5,9 m/s)
  - efficiency_vs_J.png            プロペラ効率 η vs 前進比 J
  - qprop_data/powerup_propeller_refined.prop  … 翼型確定後の QPROP 定義

使い方:
    $ python propeller_design_analysis.py
"""

import os
import sys
import glob
import subprocess
from dataclasses import dataclass

import numpy as np
import matplotlib

matplotlib.use("Agg")  # GUI 無し環境でも保存できるように
import matplotlib.pyplot as plt

# 図中の日本語が豆腐(□)にならないよう Windows 標準の日本語フォントを使う
matplotlib.rcParams["font.sans-serif"] = ["Yu Gothic", "MS Gothic", "Meiryo",
                                          "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

import propeller_config as cfg

# Windows の既定コンソール(cp932)では日本語・中点が化けるため UTF-8 に固定
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# =============================================================================
# データ構造
# =============================================================================

@dataclass
class BladeGeometry:
    """プロペラブレードの形状を保持するデータクラス。"""
    diameter: float                 # 直径 [m]
    pitch: float                    # ピッチ [m]
    n_blades: int                   # ブレード数
    radial_stations: np.ndarray     # r/R 配列(無次元)
    chord_distribution: np.ndarray  # c [m] 配列
    twist_distribution: np.ndarray  # β [rad] 配列
    airfoil_name: str               # ブレード断面翼型名

    @property
    def radius(self) -> float:
        """プロペラ半径 [m]。"""
        return self.diameter / 2.0

    @property
    def r_dim(self) -> np.ndarray:
        """各ステーションの半径 r [m]。"""
        return self.radial_stations * self.radius


# =============================================================================
# ブレード形状の構築
# =============================================================================

def load_blade_geometry_from_crazyflie(D: float = cfg.D_PROPELLER,
                                       airfoil_name: str = "S1223") -> BladeGeometry:
    """
    Crazyflie 45mm(cfnq_45p1)の c/R, β 分布を D=36mm にスケーリングして
    自作プロペラのブレード形状を構築する。

    c/R は無次元なので D を変えるだけで適用できる(設計書 v0.3 §4.3)。

    Args:
        D: 自作プロペラ直径 [m]
        airfoil_name: ブレード断面翼型名

    Returns:
        BladeGeometry インスタンス
    """
    geom = np.array(cfg.BLADE_GEOMETRY)  # (18, 3): r/R, c/R, beta[deg]
    r_R = geom[:, 0]
    c_R = geom[:, 1]
    beta_deg = geom[:, 2]

    R = D / 2.0
    chord = c_R * R                          # c [m] = (c/R) * R
    chord[0] = cfg.C_ROOT                    # 根本は強度設計の根太め化値に一元化(D5)
    twist = np.radians(beta_deg)             # β [rad]

    return BladeGeometry(
        diameter=D,
        pitch=cfg.P_PROPELLER,
        n_blades=cfg.N_BLADES,
        radial_stations=r_R,
        chord_distribution=chord,
        twist_distribution=twist,
        airfoil_name=airfoil_name,
    )


def compute_reynolds_distribution(geom: BladeGeometry, rpm: float,
                                  V: float) -> np.ndarray:
    """
    各半径での動作 Reynolds 数を計算する。

    Re = V_rel · c · ρ / μ,  V_rel = sqrt(V² + (ω·r)²)

    Args:
        geom: ブレード形状
        rpm: 回転数 [rpm]
        V: 機体速度 [m/s]

    Returns:
        各ステーションの Re 配列
    """
    omega = 2.0 * np.pi * rpm / 60.0
    v_rel = np.sqrt(V ** 2 + (omega * geom.r_dim) ** 2)
    return v_rel * geom.chord_distribution * cfg.RHO / cfg.MU


# =============================================================================
# 翼型ポーラ(XFLR5)
# =============================================================================

def _read_xflr5_polar(prefix: str, Re: int) -> tuple:
    """
    xfoil_data/propeller_blade/{prefix}_T1_Re{Re:05d}_N*.txt を 1 枚読み込む。
    NCrit 部分(N5 / N50 等)はファイルごとの表記揺れがあるためワイルドカード。

    Returns:
        (alpha, CL, CD) いずれも昇順 alpha の np.ndarray
    """
    pattern = os.path.join(cfg.PROPELLER_BLADE_POLAR_DIR,
                           f"{prefix}_T1_Re{Re:05d}_N*.txt")
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"XFLR5 ポーラ未検出: {pattern}")

    with open(matches[0], "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # 「alpha CL CD ...」ヘッダ行の 2 行先(ダッシュ行の次)から表本体
    start = None
    for i, ln in enumerate(lines):
        if ln.lstrip().lower().startswith("alpha"):
            start = i + 2
            break
    if start is None:
        raise ValueError(f"ヘッダ行が見つかりません: {matches[0]}")

    rows = []
    for ln in lines[start:]:
        parts = ln.split()
        if len(parts) < 3:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
        except ValueError:
            continue   # 非収束行・空行はスキップ
    if not rows:
        raise ValueError(f"データ行が読めません: {matches[0]}")

    rows.sort(key=lambda r: r[0])
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def load_airfoil_polars(airfoil_name: str, re_values: list,
                        alpha_deg: np.ndarray = None) -> dict:
    """
    XFLR5(v6.62)で計算済みのブレード断面翼型ポーラを読み込み、(Re, α) の
    2 次元テーブルにする。スクリーニング・QPROP 係数フィット・自作 BEM が
    共通で参照する(解答根拠を XFLR5 に一本化、§5 レビュー D4 対応)。

    各 Re のポーラを共通 α グリッドへ線形補間する。XFLR5 の収束 α 範囲
    (概ね -5〜10°)外は端点クランプ(外挿しない保守処理)。

    Args:
        airfoil_name: 翼型名(cfg.AIRFOIL_XFLR5_PREFIX のキー)
        re_values: ポーラを取得する Re のリスト
        alpha_deg: α グリッド [deg]。省略時 -5〜10° を 0.5° 刻み

    Returns:
        dict: {'name', 're_values', 'alpha', 'CL'(2D), 'CD'(2D)}
              CL[i, j] = Re=re_values[i], alpha=alpha[j] の CL
    """
    if alpha_deg is None:
        alpha_deg = np.arange(-5.0, 10.01, 0.5)

    prefix = cfg.AIRFOIL_XFLR5_PREFIX[airfoil_name]
    cl = np.zeros((len(re_values), len(alpha_deg)))
    cd = np.zeros((len(re_values), len(alpha_deg)))
    for i, re in enumerate(re_values):
        a, c_l, c_d = _read_xflr5_polar(prefix, int(re))
        cl[i, :] = np.interp(alpha_deg, a, c_l)
        cd[i, :] = np.maximum(np.interp(alpha_deg, a, c_d), 1e-4)

    return {"name": airfoil_name, "re_values": np.array(re_values, dtype=float),
            "alpha": alpha_deg, "CL": cl, "CD": cd}


def lookup_clcd(polars: dict, alpha_deg: float, re: float) -> tuple:
    """
    ポーラテーブルから (Re, α) の CL, CD を 2D 線形補間で取得する。

    範囲外の Re・α はクランプする(極低Re域の外挿を避ける保守処理)。

    Returns:
        (CL, CD)
    """
    re_v = polars["re_values"]
    a_v = polars["alpha"]
    re_c = np.clip(re, re_v[0], re_v[-1])
    a_c = np.clip(alpha_deg, a_v[0], a_v[-1])

    # Re 方向のインデックスと重み
    i = np.searchsorted(re_v, re_c) - 1
    i = np.clip(i, 0, len(re_v) - 2)
    wr = (re_c - re_v[i]) / (re_v[i + 1] - re_v[i])

    cl_lo = np.interp(a_c, a_v, polars["CL"][i])
    cl_hi = np.interp(a_c, a_v, polars["CL"][i + 1])
    cd_lo = np.interp(a_c, a_v, polars["CD"][i])
    cd_hi = np.interp(a_c, a_v, polars["CD"][i + 1])
    return (cl_lo * (1 - wr) + cl_hi * wr,
            cd_lo * (1 - wr) + cd_hi * wr)


def screen_airfoil_candidates() -> dict:
    """
    ブレード断面翼型 4 候補を動作 Re 域で NeuralFoil 解析し、L/D 最大値で
    ランキングして最良翼型を選定する(設計書 v0.3 §4.4、優先1)。

    プロペラブレードは各断面が一定の小迎角で動作するので、動作 Re 域での
    最大 L/D が高い翼型を選ぶ(効率 = 改善5項目の①)。

    Returns:
        dict: {翼型名: {'polars', 'LD_max', 'alpha_LDmax', 'CL_LDmax'}}
    """
    print("\n" + "─" * 72)
    print(" 翼型スクリーニング(ブレード断面、XFLR5 ポーラ、動作Re域)")
    print("─" * 72)
    candidates = list(cfg.AIRFOIL_XFLR5_PREFIX.keys())
    re_eval = cfg.RE_REP  # 75%R 近傍の代表 Re(設計書 §4.4)

    results = {}
    print(f"\n  {'翼型':<12}{'L/D_max':>10}{'α@LDmax':>10}{'CL@LDmax':>10}"
          f"   (代表 Re={re_eval})")
    for name in candidates:
        polars = load_airfoil_polars(name, cfg.RE_OPERATING)
        # 代表 Re での L/D を評価
        cl_row = np.array([lookup_clcd(polars, a, re_eval)[0]
                           for a in polars["alpha"]])
        cd_row = np.array([lookup_clcd(polars, a, re_eval)[1]
                           for a in polars["alpha"]])
        ld = cl_row / cd_row
        idx = int(np.argmax(ld))
        results[name] = {"polars": polars, "LD_max": ld[idx],
                         "alpha_LDmax": polars["alpha"][idx],
                         "CL_LDmax": cl_row[idx]}
        print(f"  {name:<12}{ld[idx]:>10.1f}{polars['alpha'][idx]:>10.1f}"
              f"{cl_row[idx]:>10.3f}")

    best = max(results, key=lambda k: results[k]["LD_max"])
    print(f"\n  → 最良翼型: {best}(L/D_max={results[best]['LD_max']:.1f})")
    return results, best


# =============================================================================
# 自作 BEM(ブレード要素運動量理論)
# =============================================================================

def bem_analysis(geom: BladeGeometry, rpm: float, V_inf: float,
                 polars: dict, rho: float = cfg.RHO) -> dict:
    """
    BEM(ブレード要素運動量理論)で推力・トルク・効率を計算する。

    各ブレード要素で誘導速度(軸方向 w_a・周方向 w_t)を反復解法で収束させ、
    Prandtl の翼端損失補正を加える。QPROP の交差検証用(設計書 §9.3 で
    自作 BEM は ±20% 精度想定)。

    手法:
      Va = V_inf + w_a,  Vt = ω·r - w_t
      φ  = atan2(Va, Vt),  W² = Va² + Vt²,  α = β - φ
      ブレード要素: dT = ½ρW²·B·c·(Cl·cosφ - Cd·sinφ)·dr
      運動量:       dT = 4π·r·ρ·F·Va·w_a·dr
      → w_a, w_t を反復で一致させる。

    Args:
        geom: ブレード形状
        rpm: 回転数 [rpm]
        V_inf: 機体速度 [m/s]
        polars: load_airfoil_polars() の戻り値
        rho: 空気密度 [kg/m³]

    Returns:
        dict: 'thrust'[N], 'torque'[N·m], 'power'[W], 'efficiency',
              'sectional'(各ステーションの dict)
    """
    omega = 2.0 * np.pi * rpm / 60.0
    R = geom.radius
    B = geom.n_blades
    r = geom.r_dim
    c = geom.chord_distribution
    beta = geom.twist_distribution

    # 翼端 r/R=1.0 は c→0 で寄与極小。台形積分のため全点保持
    dT = np.zeros_like(r)
    dQ = np.zeros_like(r)
    sec_alpha = np.zeros_like(r)
    sec_cl = np.zeros_like(r)
    sec_re = np.zeros_like(r)

    for k in range(len(r)):
        if c[k] < 1e-5:               # 翼端の極小弦長は寄与ゼロ扱い
            continue
        w_a, w_t = V_inf * 0.1 + 0.5, 0.5   # 初期推定
        for _ in range(200):
            Va = V_inf + w_a
            Vt = omega * r[k] - w_t
            phi = np.arctan2(Va, Vt)
            W2 = Va ** 2 + Vt ** 2
            alpha = beta[k] - phi
            Re = np.sqrt(W2) * c[k] * rho / cfg.MU
            cl, cd = lookup_clcd(polars, np.degrees(alpha), Re)
            cn = cl * np.cos(phi) - cd * np.sin(phi)   # 推力方向
            ct = cl * np.sin(phi) + cd * np.cos(phi)   # トルク方向

            # Prandtl 翼端損失
            sphi = max(abs(np.sin(phi)), 1e-3)
            f_tip = B / 2.0 * (R - r[k]) / (r[k] * sphi)
            F = max(2.0 / np.pi * np.arccos(np.exp(-min(f_tip, 30.0))), 1e-2)

            # 運動量釣り合いから誘導速度を更新(不動点反復)
            denom = 4.0 * np.pi * r[k] * rho * F * max(Va, 1e-3)
            w_a_new = 0.5 * W2 * B * c[k] * rho * cn / denom
            w_t_new = 0.5 * W2 * B * c[k] * rho * ct / denom
            # 発散防止のクランプ + 緩和
            w_a_new = np.clip(w_a_new, -0.5 * omega * R, 0.9 * omega * R)
            w_t_new = np.clip(w_t_new, -0.5 * omega * R, 0.9 * omega * R)
            if abs(w_a_new - w_a) < 1e-4 and abs(w_t_new - w_t) < 1e-4:
                w_a, w_t = w_a_new, w_t_new
                break
            w_a = 0.5 * w_a + 0.5 * w_a_new
            w_t = 0.5 * w_t + 0.5 * w_t_new

        Va = V_inf + w_a
        Vt = omega * r[k] - w_t
        phi = np.arctan2(Va, Vt)
        W2 = Va ** 2 + Vt ** 2
        alpha = beta[k] - phi
        Re = np.sqrt(W2) * c[k] * rho / cfg.MU
        cl, cd = lookup_clcd(polars, np.degrees(alpha), Re)
        cn = cl * np.cos(phi) - cd * np.sin(phi)
        ct = cl * np.sin(phi) + cd * np.cos(phi)
        dT[k] = 0.5 * rho * W2 * B * c[k] * cn
        dQ[k] = 0.5 * rho * W2 * B * c[k] * ct * r[k]
        sec_alpha[k] = np.degrees(alpha)
        sec_cl[k] = cl
        sec_re[k] = Re

    thrust = float(np.trapezoid(dT, r))
    torque = float(np.trapezoid(dQ, r))
    power = torque * omega
    eff = (thrust * V_inf / power) if power > 1e-9 and V_inf > 0 else 0.0

    return {"thrust": thrust, "torque": torque, "power": power,
            "efficiency": eff, "rpm": rpm, "V": V_inf,
            "sectional": {"r_R": geom.radial_stations, "dT": dT, "dQ": dQ,
                          "alpha": sec_alpha, "CL": sec_cl, "Re": sec_re}}


# =============================================================================
# QPROP 連携(主たる数値、handoff §6.2)
# =============================================================================

def run_qprop(prop_file: str, motor_file: str, V: float, rpm: float,
              volt: float = 0.0) -> dict:
    """
    QPROP 1.22 を subprocess で呼び、動作点の性能を取得する。

    RPM を直接指定するため推力・トルクはモーターパラメータに非依存
    (handoff §2.3 注記)。qprop.exe は静的リンク済みで PATH 非依存。

    重要(Windows 移植時のハマりどころ):
      QPROP(Fortran)のファイル名バッファは 48 文字しかなく、長い絶対パスは
      途中で切られて「Prop file not found → 既定プロペラ」になる。これを避け、
      cwd を qprop_data/ に固定してファイル名のみを渡す。

    Args:
        prop_file: .prop ファイルパス(絶対/相対どちらでも可、basename を使用)
        motor_file: .mot ファイルパス(同上)
        V: 機体速度 [m/s]
        rpm: 回転数 [rpm]
        volt: 電圧 [V](0 で RPM 指定モード)

    Returns:
        dict: 'thrust'[N], 'torque'[N·m], 'power'[W], 'eff_prop',
              'rpm'[rpm], 'CT', 'CP', 'V'[m/s]
    """
    result = subprocess.run(
        [cfg.QPROP_BIN, os.path.basename(prop_file),
         os.path.basename(motor_file), str(V), str(rpm), str(volt)],
        capture_output=True, text=True, timeout=60, cwd=cfg.QPROP_DATA_DIR,
    )
    if result.returncode != 0:
        raise RuntimeError(f"QPROP 異常終了:\n{result.stderr}\n{result.stdout}")
    if "not found" in result.stdout:
        raise RuntimeError(f"QPROP がファイルを読めません:\n{result.stdout[:300]}")

    # 出力をパース: "V(m/s) rpm Dbeta T(N) ..." ヘッダ直後の数値行
    lines = result.stdout.splitlines()
    header_idx = None
    for i, ln in enumerate(lines):
        if "V(m/s)" in ln and "rpm" in ln:
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError(f"QPROP 出力をパースできません:\n{result.stdout}")

    # ヘッダ直後の最初の数値データ行(行頭 '#' は除去)
    for ln in lines[header_idx + 1:]:
        toks = ln.lstrip("#").split()
        if len(toks) >= 17:
            try:
                vals = [float(t) for t in toks[:17]]
            except ValueError:
                continue
            return {"V": vals[0], "rpm": vals[1], "thrust": vals[3],
                    "torque": vals[4], "power": vals[5],
                    "eff_prop": vals[9], "CT": vals[11], "CP": vals[12]}
    raise RuntimeError(f"QPROP データ行が見つかりません:\n{result.stdout}")


def qprop_sweep_rpm(prop_file: str, motor_file: str, V: float,
                    rpm_list) -> dict:
    """指定速度 V で RPM を振り、推力・効率を取得する。"""
    T, eff = [], []
    for rpm in rpm_list:
        res = run_qprop(prop_file, motor_file, V, rpm)
        T.append(res["thrust"])
        eff.append(res["eff_prop"])
    return {"rpm": np.array(rpm_list, dtype=float),
            "thrust": np.array(T), "eff": np.array(eff)}


# =============================================================================
# QPROP 翼型パラメータの精緻化(優先1: .prop 更新)
# =============================================================================

def fit_qprop_airfoil_params(polars: dict) -> dict:
    """
    NeuralFoil ポーラから QPROP の単一翼型モデル係数をフィットする。

    QPROP 翼型モデル:
      Cl = CL0 + CL_a·α        (CLmin..CLmax でクランプ)
      Cd = CD0 + CD2·(Cl-CLCD0)²   (Cl>CLCD0 は CD2u、Cl<CLCD0 は CD2l)
      Re 補正: Cd ∝ (Re/REref)^REexp

    動作 Re 域の代表(16,000)でフィットする。

    Returns:
        dict: QPROP .prop に書く翼型係数
    """
    re_ref = float(cfg.RE_REP)
    a_deg = polars["alpha"]
    cl = np.array([lookup_clcd(polars, a, re_ref)[0] for a in a_deg])
    cd = np.array([lookup_clcd(polars, a, re_ref)[1] for a in a_deg])
    a_rad = np.radians(a_deg)

    # 線形域(-3〜+6°)で CL0, CL_a を最小二乗フィット
    lin = (a_deg >= -3.0) & (a_deg <= 6.0)
    slope, intercept = np.polyfit(a_rad[lin], cl[lin], 1)  # [1/rad]
    cl0, cl_a = float(intercept), float(slope)

    cl_min, cl_max = float(cl.min()), float(cl.max())

    # ドラッグポーラ Cd = CD0 + CD2·(Cl-CLCD0)²
    i_min = int(np.argmin(cd))
    cd0 = float(cd[i_min])
    cl_cd0 = float(cl[i_min])
    upper = cl > cl_cd0
    lower = cl < cl_cd0
    cd2u = float(np.polyfit(cl[upper] - cl_cd0, cd[upper] - cd0, 2)[0]) \
        if upper.sum() >= 3 else 0.04
    cd2l = float(np.polyfit(cl[lower] - cl_cd0, cd[lower] - cd0, 2)[0]) \
        if lower.sum() >= 3 else 0.04
    cd2u = float(np.clip(cd2u, 0.005, 0.5))
    cd2l = float(np.clip(cd2l, 0.005, 0.5))

    return {"CL0": cl0, "CL_a": cl_a, "CLmin": cl_min, "CLmax": cl_max,
            "CD0": cd0, "CD2u": cd2u, "CD2l": cd2l, "CLCD0": cl_cd0,
            "REref": re_ref, "REexp": -0.5}


def write_refined_prop(geom: BladeGeometry, params: dict,
                       out_path: str) -> str:
    """
    翼型確定後の QPROP .prop ファイルを書き出す(優先1: .prop 更新)。

    形状(18 ステーション)は凍結値、翼型係数のみ NeuralFoil フィット値に
    更新する。

    Returns:
        書き出したファイルパス
    """
    lines = []
    lines.append("")
    lines.append(f"POWERUP 4.0 self-made propeller D=36mm P=27mm "
                 f"({geom.airfoil_name} blade section, XFLR5-refined)")
    lines.append("")
    lines.append(f" {geom.n_blades}     {geom.radius:.4f}   "
                 f"! Nblades   R [m]")
    lines.append("")
    lines.append(f"! Blade airfoil: {geom.airfoil_name}, "
                 f"XFLR5 polar fit @ Re={params['REref']:.0f}")
    lines.append(f" {params['CL0']:.4f}   {params['CL_a']:.4f}   "
                 f"! CL0       CL_a  (per rad)")
    lines.append(f" {params['CLmin']:.4f}   {params['CLmax']:.4f}   "
                 f"! CLmin     CLmax")
    lines.append("")
    lines.append(f" {params['CD0']:.5f}  {params['CD2u']:.5f}  "
                 f"{params['CD2l']:.5f}  {params['CLCD0']:.4f}   "
                 f"! CD0   CD2u   CD2l   CLCD0")
    lines.append(f" {params['REref']:.0f}  {params['REexp']:.2f}               "
                 f"! REref REexp")
    lines.append("")
    lines.append(" 0.018   0.001   1.0    ! Rfac    Cfac    Bfac")
    lines.append(" 0.0     0.0     0.0    ! Radd    Cadd    Badd")
    lines.append("")
    lines.append("#  r/R    chord[mm]   beta[deg]")
    for r_R, c, b in zip(geom.radial_stations,
                         geom.chord_distribution,
                         geom.twist_distribution):
        lines.append(f" {r_R:.2f}    {c * 1000:.2f}         "
                     f"{np.degrees(b):.1f}")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out_path


# =============================================================================
# プロット
# =============================================================================

def plot_airfoil_polars(screen_results: dict, best: str, save_path: str):
    """翼型 4 候補のポーラ比較(計算書図 #4)。代表 Re=15,000。"""
    cfg.ensure_plots_dir()
    re_eval = cfg.RE_REP
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    colors = {"S1223": "#d32f2f", "E205": "#1976d2", "AG12": "#388e3c",
              "camber4_thick2": "#f57c00"}
    for name, res in screen_results.items():
        p = res["polars"]
        cl = np.array([lookup_clcd(p, a, re_eval)[0] for a in p["alpha"]])
        cd = np.array([lookup_clcd(p, a, re_eval)[1] for a in p["alpha"]])
        lw = 2.4 if name == best else 1.3
        col = colors.get(name, "gray")
        lab = f"{name}（採用）" if name == best else name
        axes[0].plot(p["alpha"], cl, color=col, lw=lw, label=lab)
        axes[1].plot(p["alpha"], cd, color=col, lw=lw, label=lab)
        axes[2].plot(p["alpha"], cl / cd, color=col, lw=lw, label=lab)
    axes[0].set(xlabel="α [deg]", ylabel="CL", title="CL vs α")
    axes[1].set(xlabel="α [deg]", ylabel="CD", title="CD vs α")
    axes[2].set(xlabel="α [deg]", ylabel="L/D", title="L/D vs α")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"Blade-section airfoil candidates  (XFLR5, Re={re_eval})",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


def plot_blade_distributions(geom: BladeGeometry, save_chord: str,
                             save_twist: str):
    """弦長分布・ねじり角分布(自作 vs Crazyflie 原寸)。計算書図 #2, #3。"""
    cfg.ensure_plots_dir()
    cf = np.array(cfg.BLADE_GEOMETRY)
    r_R = cf[:, 0]
    cf_chord_mm = cf[:, 1] * (0.045 / 2) * 1000   # Crazyflie 原寸 D=45mm
    self_chord_mm = geom.chord_distribution * 1000

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(r_R, self_chord_mm, "o-", color="#d32f2f", lw=2,
            label="自作 D=36mm")
    ax.plot(r_R, cf_chord_mm, "s--", color="#1976d2", lw=1.5,
            label="Crazyflie D=45mm（原寸）")
    ax.set(xlabel="r/R", ylabel="chord c [mm]",
           title="Blade chord distribution")
    ax.grid(alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_chord, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_chord}")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(r_R, np.degrees(geom.twist_distribution), "o-",
            color="#d32f2f", lw=2, label="自作 D=36mm")
    ax.plot(r_R, cf[:, 2], "s--", color="#1976d2", lw=1.5,
            label="Crazyflie cfnq_45p1")
    ax.set(xlabel="r/R", ylabel="twist β [deg]",
           title="Blade twist distribution")
    ax.grid(alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_twist, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_twist}")


def plot_geometry_3view(geom: BladeGeometry, save_path: str):
    """プロペラ三面図(正面・側面・上面)。計算書図 #1。"""
    cfg.ensure_plots_dir()
    r = geom.r_dim * 1000        # mm
    c = geom.chord_distribution * 1000
    R = geom.radius * 1000

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # --- 上面図(平面形、2 枚ブレード)---
    ax = axes[0]
    for sign in (1, -1):
        ax.fill_between(sign * r, -c / 2, c / 2, color="#bcd4e6",
                        edgecolor="black", lw=0.8)
    hub = plt.Circle((0, 0), cfg.HUB_OUTER_DIAMETER * 1000 / 2,
                     color="gray", ec="black")
    ax.add_patch(hub)
    ax.set(xlabel="r [mm]", ylabel="chord [mm]", title="上面図(平面形)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)

    # --- 正面図(回転投影円)---
    ax = axes[1]
    circle = plt.Circle((0, 0), R, fill=False, ec="#1976d2", lw=2)
    ax.add_patch(circle)
    ax.add_patch(plt.Circle((0, 0), cfg.HUB_OUTER_DIAMETER * 1000 / 2,
                            color="gray", ec="black"))
    ax.plot([-R, R], [0, 0], "k--", lw=0.8)
    ax.plot([0, 0], [-R, R], "k--", lw=0.8)
    ax.text(0, R * 0.6, f"D = {geom.diameter * 1000:.0f} mm",
            ha="center", fontsize=10)
    ax.set(xlim=(-R * 1.2, R * 1.2), ylim=(-R * 1.2, R * 1.2),
           title="正面図(回転投影)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)

    # --- 側面図(ねじり角)---
    ax = axes[2]
    beta_deg = np.degrees(geom.twist_distribution)
    ax.plot(geom.radial_stations, beta_deg, "o-", color="#d32f2f", lw=2)
    ax.set(xlabel="r/R", ylabel="twist β [deg]",
           title=f"側面図相当: ねじり角  (P={geom.pitch * 1000:.0f}mm)")
    ax.grid(alpha=0.3)

    fig.suptitle(f"Self-made propeller 3-view  "
                 f"(D={geom.diameter * 1000:.0f}mm, P={geom.pitch * 1000:.0f}mm, "
                 f"{geom.n_blades} blades, {geom.airfoil_name})",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


def plot_thrust_vs_rpm(sweeps: dict, save_path: str):
    """推力 vs RPM 曲線(V=0,5,7.5,9 m/s)。計算書図 #5。"""
    cfg.ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(7.5, 5))
    colors = {0.0: "#888", 5.0: "#388e3c", 7.5: "#1976d2", 9.0: "#d32f2f"}
    for V, sw in sweeps.items():
        ax.plot(sw["rpm"] / 1000, sw["thrust"] * 1000, "o-",
                color=colors.get(V, "k"), label=f"V = {V} m/s")
    ax.axhline(cfg.T_REQUIRED_PER_MOTOR * 1000, color="red", ls=":",
               label=f"必要推力 {cfg.T_REQUIRED_PER_MOTOR * 1000:.1f} mN")
    ax.axvline(cfg.RPM_CRUISE / 1000, color="gray", ls="--", alpha=0.6)
    ax.set(xlabel="RPM [×1000]", ylabel="thrust T [mN]",
           title="Thrust vs RPM  (QPROP)")
    ax.grid(alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


def plot_efficiency_vs_J(J: np.ndarray, eff: np.ndarray, J_cruise: float,
                         save_path: str):
    """プロペラ効率 η vs 前進比 J 曲線。計算書図 #6。"""
    cfg.ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(J, eff, "o-", color="#1976d2", lw=2)
    ax.axvline(J_cruise, color="red", ls="--",
               label=f"巡航 J = {J_cruise:.3f}")
    i_best = int(np.argmax(eff))
    ax.plot(J[i_best], eff[i_best], "r*", ms=16,
            label=f"η_max = {eff[i_best]:.3f} @ J={J[i_best]:.3f}")
    ax.set(xlabel="advance ratio  J = V/(nD)", ylabel="propeller efficiency η",
           title="Efficiency vs advance ratio  (QPROP, RPM=22,100)")
    ax.grid(alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


# =============================================================================
# メイン
# =============================================================================

def main():
    print("=" * 72)
    print(" §5 プロペラ設計  空力性能解析  (作成 2026-05-16)")
    print("=" * 72)
    print(f"\n[設計諸元]")
    print(f"  直径 D            = {cfg.D_PROPELLER * 1000:.0f} mm")
    print(f"  ピッチ P          = {cfg.P_PROPELLER * 1000:.0f} mm  "
          f"(P/D = {cfg.P_OVER_D:.2f})")
    print(f"  ブレード数        = {cfg.N_BLADES}")
    print(f"  巡航動作 RPM      = {cfg.RPM_CRUISE} rpm")
    print(f"  1モーター必要推力 = {cfg.T_REQUIRED_PER_MOTOR * 1000:.2f} mN")

    # ---------- 翼型スクリーニング(優先1)----------
    screen_results, best_airfoil = screen_airfoil_candidates()

    # ---------- ブレード形状 ----------
    geom = load_blade_geometry_from_crazyflie(D=cfg.D_PROPELLER,
                                              airfoil_name=best_airfoil)
    print(f"\n[ブレード形状]  Crazyflie cfnq_45p1 を D=36mm にスケーリング")
    print(f"  最大弦長 = {geom.chord_distribution.max() * 1000:.2f} mm "
          f"(r/R={geom.radial_stations[np.argmax(geom.chord_distribution)]:.2f})")
    print(f"  ねじり角 = 根本 {np.degrees(geom.twist_distribution[0]):.1f}° "
          f"→ 先端 {np.degrees(geom.twist_distribution[-1]):.1f}°")

    # ---------- 動作 Re 域 ----------
    re_dist = compute_reynolds_distribution(geom, cfg.RPM_CRUISE, cfg.V_CRUISE)
    # 翼端 c→0 を除いた有効域
    valid = geom.chord_distribution > 1e-4
    print(f"\n[動作 Re 域 @ {cfg.RPM_CRUISE}rpm, V={cfg.V_CRUISE}m/s]")
    print(f"  Re = {re_dist[valid].min():.0f} 〜 {re_dist[valid].max():.0f}")

    # ---------- QPROP 翼型パラメータ精緻化 → .prop 更新(優先1)----------
    print("\n" + "─" * 72)
    print(" QPROP 翼型パラメータの精緻化(.prop 更新)")
    print("─" * 72)
    params = fit_qprop_airfoil_params(screen_results[best_airfoil]["polars"])
    print(f"  翼型 {best_airfoil} の XFLR5 フィット結果(Re={cfg.RE_REP:,}):")
    print(f"    CL0={params['CL0']:.3f}  CL_a={params['CL_a']:.2f}/rad  "
          f"CLmax={params['CLmax']:.2f}")
    print(f"    CD0={params['CD0']:.4f}  CD2u={params['CD2u']:.3f}  "
          f"CD2l={params['CD2l']:.3f}  CLCD0={params['CLCD0']:.3f}")
    refined_prop = os.path.join(cfg.QPROP_DATA_DIR,
                                "powerup_propeller_refined.prop")
    write_refined_prop(geom, params, refined_prop)
    print(f"  → 保存: {refined_prop}")

    motor_file = os.path.join(cfg.QPROP_DATA_DIR, "powerup_motor.mot")

    # ---------- 設計点の QPROP vs 自作 BEM 交差検証 ----------
    print("\n" + "─" * 72)
    print(" 設計点性能(巡航 22,100rpm / V=7.5m/s):QPROP vs 自作BEM")
    print("─" * 72)
    qp = run_qprop(refined_prop, motor_file, cfg.V_CRUISE, cfg.RPM_CRUISE)
    bem = bem_analysis(geom, cfg.RPM_CRUISE, cfg.V_CRUISE,
                       screen_results[best_airfoil]["polars"])
    print(f"\n  {'量':<22}{'QPROP':>14}{'自作BEM':>14}")
    print(f"  {'推力 T [mN]':<22}{qp['thrust'] * 1000:>14.2f}"
          f"{bem['thrust'] * 1000:>14.2f}")
    print(f"  {'トルク Q [mN·m]':<20}{qp['torque'] * 1000:>14.3f}"
          f"{bem['torque'] * 1000:>14.3f}")
    print(f"  {'プロペラ効率 η':<21}{qp['eff_prop']:>14.3f}"
          f"{bem['efficiency']:>14.3f}")
    diff = abs(qp["thrust"] - bem["thrust"]) / qp["thrust"] * 100
    print(f"\n  推力の QPROP-BEM 乖離 = {diff:.1f} %  "
          f"({'✓ ±20%以内' if diff <= 20 else '△ 要確認'})")
    print(f"  必要推力 26.5 mN に対し QPROP 予測は "
          f"{qp['thrust'] * 1000 / (cfg.T_REQUIRED_PER_MOTOR * 1000) * 100:.0f}%")

    # ---------- 推力 vs RPM 曲線(QPROP)----------
    print("\n" + "─" * 72)
    print(" 推力曲線・効率曲線の生成(QPROP スイープ)")
    print("─" * 72)
    rpm_list = [15000, 20000, 22100, 26600, 30000, 35000, 40000]
    sweeps = {}
    for V in (0.0, 5.0, 7.5, 9.0):
        sweeps[V] = qprop_sweep_rpm(refined_prop, motor_file, V, rpm_list)
    t_static_max = sweeps[0.0]["thrust"][-1]   # V=0, 40K
    print(f"  最大静推力(40K rpm, V=0)= {t_static_max * 1000:.1f} mN")

    # 巡航必要推力に到達する RPM(V=7.5)を細かいスイープで特定
    fine_rpm = np.arange(20000.0, 34001.0, 1000.0)
    fine = qprop_sweep_rpm(refined_prop, motor_file, cfg.V_CRUISE, fine_rpm)
    t_req = cfg.T_REQUIRED_PER_MOTOR
    if fine["thrust"].min() <= t_req <= fine["thrust"].max():
        rpm_cruise_thrust = float(np.interp(t_req, fine["thrust"], fine["rpm"]))
        print(f"  巡航必要推力 {t_req * 1000:.1f} mN 到達 RPM(V=7.5)"
              f"= {rpm_cruise_thrust:.0f} rpm  "
              f"(モーター上限 {cfg.RPM_MAX} rpm の "
              f"{rpm_cruise_thrust / cfg.RPM_MAX * 100:.0f}%)")
        qp_cruise = run_qprop(refined_prop, motor_file, cfg.V_CRUISE,
                              rpm_cruise_thrust)
        print(f"  巡航動作点({rpm_cruise_thrust:.0f}rpm, V=7.5)"
              f": T={qp_cruise['thrust'] * 1000:.1f} mN, "
              f"η={qp_cruise['eff_prop']:.3f}")
    else:
        rpm_cruise_thrust = None
        print(f"  [注意] 必要推力 {t_req * 1000:.1f} mN は "
              f"{fine_rpm[0]:.0f}-{fine_rpm[-1]:.0f} rpm の範囲外")

    # ---------- 効率 vs 前進比 J ----------
    # 高前進比(推力≈0 近傍)では QPROP の effprop が破綻するため、
    # 推力 > 0 かつ η が物理範囲(0〜1)の点のみ採用する。
    V_range = np.linspace(0.5, 11.0, 22)
    n = cfg.RPM_CRUISE / 60.0
    J = V_range / (n * cfg.D_PROPELLER)
    eff_J, thr_J = [], []
    for V in V_range:
        res = run_qprop(refined_prop, motor_file, V, cfg.RPM_CRUISE)
        eff_J.append(res["eff_prop"])
        thr_J.append(res["thrust"])
    eff_J = np.array(eff_J)
    thr_J = np.array(thr_J)
    physical = (thr_J > 0) & (eff_J > 0) & (eff_J < 1.0)
    J_cruise = cfg.V_CRUISE / (n * cfg.D_PROPELLER)
    i_best = int(np.where(physical)[0][np.argmax(eff_J[physical])])
    print(f"  巡航前進比 J = {J_cruise:.3f}")
    print(f"  最大効率 η = {eff_J[i_best]:.3f} @ J = {J[i_best]:.3f}")
    print(f"  推力ゼロ到達 J ≈ {J[thr_J > 0][-1]:.3f}（これ以上は windmill 域)")

    # ---------- プロット生成 ----------
    print("\n" + "─" * 72)
    print(" 図表生成(plots/)")
    print("─" * 72)
    plot_geometry_3view(geom, cfg.plot_path("propeller_3view.png"))
    plot_blade_distributions(geom,
                             cfg.plot_path("blade_chord_distribution.png"),
                             cfg.plot_path("blade_twist_distribution.png"))
    plot_airfoil_polars(screen_results, best_airfoil,
                        cfg.plot_path("airfoil_polars.png"))
    plot_thrust_vs_rpm(sweeps, cfg.plot_path("thrust_vs_rpm_curves.png"))
    plot_efficiency_vs_J(J[physical], eff_J[physical], J_cruise,
                         cfg.plot_path("efficiency_vs_J.png"))

    print("\n" + "=" * 72)
    print(" 確定値サマリー(計算書 §5 提供物)")
    print("=" * 72)
    print(f"  ブレード断面翼型      = {best_airfoil}")
    print(f"  巡航推力(QPROP)     = {qp['thrust'] * 1000:.2f} mN "
          f"(必要 26.5 mN)")
    print(f"  巡航効率(QPROP)     = {qp['eff_prop']:.3f}")
    print(f"  最大効率              = {eff_J[i_best]:.3f} @ J={J[i_best]:.3f}")
    print(f"  最大静推力(40K rpm) = {t_static_max * 1000:.1f} mN")
    if rpm_cruise_thrust:
        print(f"  巡航必要推力到達RPM   = {rpm_cruise_thrust:.0f} rpm (V=7.5)")
    print("\n[完了]")

    return {"geom": geom, "best_airfoil": best_airfoil, "qprop_design": qp,
            "bem_design": bem, "sweeps": sweeps, "eff_J": eff_J, "J": J,
            "refined_prop": refined_prop}


if __name__ == "__main__":
    main()
