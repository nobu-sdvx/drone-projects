"""
2026-05-17_boom_attachment_analysis.py
=======================================

ツインブーム接合方式変更による影響解析(数値計算のみ)

2026-05-17 設計レビューで決定された新方式の重量・重心・強度への影響を一括算出する。

  旧方式: ブームが主翼を弦方向に貫通(全長 約 200 mm)
  新方式: ブームを主翼後縁から 30 mm 深埋め込み + EPS スリーブで根元補強

本スクリプトは数値計算のみを行い、計算書(§2 主翼・§3 尾翼・§7 重量重心)の
文書修正は行わない(後続ハンドオフで実施)。

実行:  python 2026-05-17_boom_attachment_analysis.py
出力:  標準出力に Task 1〜7 の式・入力値・結果

注記: 本解析に数値積分は無いため numpy は使用しない(規約上 np.trapz は使わず、
       積分が必要な場合のみ np.trapezoid を使う方針)。
"""

import math

# =============================================================================
# 定数定義(凍結値、変更不可)
# =============================================================================

# --- 機体ジオメトリ(重量重心書 v3 §7、主翼書 §2)---
WING_LE = 128.0          # 主翼前縁位置(設計点)[mm]  出典: 重量重心書 v3 §7.5
WING_TE = 188.0          # 主翼後縁位置(設計点)[mm]  = WING_LE + 60
C_WING = 60.0            # 主翼弦 [mm]                出典: §2
X_WING_SPAR = 146.0      # 主翼スパー位置 [mm]        弦の30%(前縁から18mm)、AG36 最大厚位置
X_TAIL_END = 305.0       # 尾翼後端位置 [mm]          出典: §7
M_TOTAL_OLD = 36.47      # 機体総重量(現行)[g]      出典: §7.2 詳細表
X_CG_OLD = 148.94        # 設計重心 x_cg(現行)[mm]  出典: §7.4
SUM_MX_OLD = 5431.8      # Σ(m·x)(現行)[g·mm]       出典: §7.2 詳細表

# --- ブーム新仕様 ---
BOOM_LINEAR_DENSITY = 4.87e-3   # 線密度 [g/mm]  ρ1.55g/cm³ × A=π·(1mm)²
BOOM_DIAMETER = 2.0             # φ2 mm CFRP ロッド [mm]
N_BOOM = 2                      # 本数(左右対称)
BOOM_EMBED = 30.0               # 主翼内 埋め込み長 [mm]
X_BOOM_FRONT = WING_TE - BOOM_EMBED   # ブーム前端 = 188 − 30 = 158 mm
X_BOOM_REAR = X_TAIL_END              # ブーム後端 = 305 mm
BOOM_LENGTH = X_BOOM_REAR - X_BOOM_FRONT          # 全長 = 147 mm
BOOM_CANTILEVER = X_TAIL_END - WING_TE            # 空中露出(片持ち梁)長 = 117 mm
X_BOOM_CG = (X_BOOM_FRONT + X_BOOM_REAR) / 2.0    # ブーム CG = 231.5 mm(線密度一定)
M_BOOM_OLD = 1.8                # 旧ブーム合計質量 [g]  出典: §7.2 現行値
X_BOOM_OLD = 213.0              # 旧ブーム CG 位置 [mm]  出典: §7.2 現行値

L_WORST = 0.100                 # 設計荷重(最悪)[N]  出典: 尾翼書 §3.1.2
L_PER_BOOM = L_WORST / 2.0      # 1本あたり負担 [N]   左右均等仮定
SIGMA_ULT_CFRP = 1500e6         # CFRP 引張強度 [Pa]   出典: 主翼書 §2.4

# --- EPS スリーブ新仕様(新規部品)---
SLEEVE_OD = 5.0          # 外径 [mm]
SLEEVE_ID = 2.1          # 内径 [mm](ブーム嵌合公差込み、主翼貫通穴と同径)
SLEEVE_LENGTH = 15.0     # 長さ [mm](主翼後縁から後方への突出区間)
SLEEVE_DENSITY = 0.025e-3  # EPS 密度 [g/mm³]  = 0.025 g/cm³
N_SLEEVE = 2             # 個数(左右対称)
X_SLEEVE_CG = (WING_TE + (WING_TE + SLEEVE_LENGTH)) / 2.0   # = 195.5 mm

# --- 接着剤(エポキシ系)---
TAU_EPOXY = 10e6         # せん断強度(保守値)[Pa]  出典: プロペラ書 §5.5
ETA_ADHESIVE = 0.10      # 有効率                     出典: プロペラ書 §5.5
D_HOLE = 2.1             # 主翼貫通穴径 [mm]          スリーブ ID と同径
F_AXIAL = 0.05           # 軸方向引き抜き力(暫定値)[N]

# --- EPS 物性 ---
SIGMA_EPS_Y = 0.4e6      # EPS 圧縮降伏応力 [Pa]      出典: 尾翼書 §4.3(20倍EPS保守値)
E_EPS = 5e6              # EPS 弾性率 [Pa]            出典: 尾翼書 §4.3

# --- 判定目標 ---
SF_TARGET = 4.0          # 強度安全率 目標
SM_RANGE = (5.0, 15.0)   # 静安定余裕 目標範囲 [%]
VH_RANGE = (0.30, 0.60)  # 水平尾翼容積比 目標範囲

# --- 静安定の整合済み v3 入力値(整合済み概略計算書 §3/§7、2026-05-17)---
AR_WING = 6.0            # 主翼アスペクト比
AR_H = 148.0 / 36.0      # 水平尾翼アスペクト比 = 4.111(148×36mm、§3 U字版v3)
VH_TARGET = 0.55         # 水平尾翼容積比 設計目標
SH_CM2 = 53.28           # 水平尾翼面積 [cm²](§3 U字版v3)
ELL_H = 134.25           # 主翼AC〜尾翼AC 距離 [mm](§3/§7、x_cg に非依存)
S_WING_CM2 = 216.0       # 主翼面積 [cm²]
SM_BASELINE = 14.7       # 整合済み §7.6.2 設計点 SM [%](重心を翼弦35%に置いた定義値)


# =============================================================================
# 関数群
# =============================================================================

def compute_boom_mass():
    """
    Task 1 — ブーム新質量。

    入力: 線密度 [g/mm]、全長 [mm](定数)
    出力: dict(1本質量、合計質量、旧値からの差分)[g]
    前提: 線密度一定(φ2mm CFRP ロッド、ρ=1.55 g/cm³)
    """
    m_per_boom = BOOM_LINEAR_DENSITY * BOOM_LENGTH
    m_total = m_per_boom * N_BOOM
    delta = m_total - M_BOOM_OLD
    return {"per_boom": m_per_boom, "total": m_total, "delta": delta}


def compute_sleeve_mass():
    """
    Task 2 — EPS スリーブ質量。

    入力: 外径・内径・長さ・密度(定数)
    出力: dict(実効体積、1個質量、合計質量)
    前提: 中空円筒(外径円筒 − 内径円筒)、密度 0.025 g/cm³
    """
    v_outer = math.pi * (SLEEVE_OD / 2.0) ** 2 * SLEEVE_LENGTH
    v_inner = math.pi * (SLEEVE_ID / 2.0) ** 2 * SLEEVE_LENGTH
    v_sleeve = v_outer - v_inner
    m_per = v_sleeve * SLEEVE_DENSITY
    m_total = m_per * N_SLEEVE
    return {"v_outer": v_outer, "v_inner": v_inner, "v_sleeve": v_sleeve,
            "per_sleeve": m_per, "total": m_total}


def section_modulus_circular(d_mm):
    """円形断面の断面係数 W = π·d³/32 [mm³]。"""
    return math.pi * d_mm ** 3 / 32.0


def evaluate_boom_strength():
    """
    Task 3 — ブーム空中露出区間の強度(片持ち梁モデル)。

    入力: 1本負担荷重 [N]、片持ち長 [m]、断面係数(定数)
    出力: dict(曲げモーメント、応力、安全率、旧SFとの比較)
    前提: 空中露出長 117 mm を片持ち梁長とする。埋め込み区間は梁に含めない。
    """
    ell_m = BOOM_CANTILEVER / 1000.0
    m_bend = L_PER_BOOM * ell_m                       # [N·m]
    w_s_m3 = section_modulus_circular(BOOM_DIAMETER) * 1e-9   # [m³]
    sigma = m_bend / w_s_m3                           # [Pa]
    sf = SIGMA_ULT_CFRP / sigma
    return {"M": m_bend, "W_s": w_s_m3, "sigma": sigma, "SF": sf,
            "SF_old": 118.0}


def evaluate_pullout():
    """
    Task 4a — 埋め込み部の軸方向引き抜き保持力。

    入力: 穴径・埋め込み長・接着強度・有効率・軸力(定数)
    出力: dict(接着面積、有効せん断強度、保持力、安全率)
    前提: 接着面は φ2.1mm 穴 × 埋め込み長 30mm の円筒側面。
    """
    a_adhesive = math.pi * D_HOLE * BOOM_EMBED          # [mm²]
    a_adhesive_m2 = a_adhesive * 1e-6                   # [m²]
    tau_eff = TAU_EPOXY * ETA_ADHESIVE                  # [Pa]
    f_pullout = tau_eff * a_adhesive_m2                 # [N]
    sf = f_pullout / F_AXIAL
    return {"A": a_adhesive, "tau_eff": tau_eff, "F_pullout": f_pullout,
            "F_axial": F_AXIAL, "SF": sf}


def evaluate_rotation_mode(m_embed):
    """
    Task 4b — 埋め込み部の曲げによる引き抜き(回転モード)。

    入力: 埋め込み部に作用する曲げモーメント m_embed [N·m]
    出力: dict(反力作用幅、偶力、EPS圧縮応力、安全率)
    前提: 反力作用幅 a = (2/3)·埋め込み長(保守的仮定)。
          偶力 F = M/a が EPS を a×穴径 の面で圧縮する。
    """
    a_react = (2.0 / 3.0) * BOOM_EMBED                  # [mm]
    a_react_m = a_react / 1000.0
    f_couple = m_embed / a_react_m                      # [N]
    area_m2 = a_react_m * (D_HOLE / 1000.0)             # [m²]
    sigma_eps = f_couple / area_m2                      # [Pa]
    sf = SIGMA_EPS_Y / sigma_eps
    return {"a": a_react, "F_couple": f_couple, "sigma_eps": sigma_eps, "SF": sf}


def recompute_cg(boom, sleeve):
    """
    Task 5 — 重心位置の再計算。

    入力: ブーム質量 dict、スリーブ質量 dict
    出力: dict(新総質量、新Σmx、新重心)
    前提: ブーム以外・スリーブ以外の部品は §7.2 から不変。
          Σmx の更新 = 旧Σmx − 旧ブームモーメント + 新ブームモーメント + スリーブモーメント。
    """
    m_total_new = M_TOTAL_OLD + boom["delta"] + sleeve["total"]
    sum_mx_new = (SUM_MX_OLD
                  - M_BOOM_OLD * X_BOOM_OLD
                  + boom["total"] * X_BOOM_CG
                  + sleeve["total"] * X_SLEEVE_CG)
    x_cg_new = sum_mx_new / m_total_new
    return {"m_total": m_total_new, "sum_mx": sum_mx_new, "x_cg": x_cg_new,
            "delta_x_cg": x_cg_new - X_CG_OLD}


def neutral_point_frac():
    """中立点位置 x_NP/c(Drela MIT 16.01 "Lab 8 Notes" 近似式)。"""
    return (0.25 + (1.0 - 4.0 / (AR_WING + 2.0)) * VH_TARGET
            * (1.0 + 2.0 / AR_WING) / (1.0 + 2.0 / AR_H))


def check_stability(cg):
    """
    Task 6 — 静安定性指標の確認。

    入力: 重心 dict(recompute_cg の出力)
    出力: dict(中立点、SM、Vh と各判定)
    前提: 主翼配置 wing_LE=128.0 据え置き。中立点式・Vh は x_cg に非依存。
    """
    xnp_frac = neutral_point_frac()
    x_np = WING_LE + xnp_frac * C_WING
    sm = (x_np - cg["x_cg"]) / C_WING * 100.0           # [%]
    vh = (SH_CM2 * ELL_H) / (S_WING_CM2 * C_WING)
    return {"xnp_frac": xnp_frac, "x_np": x_np, "SM": sm, "Vh": vh,
            "SM_ok": SM_RANGE[0] <= sm <= SM_RANGE[1],
            "Vh_ok": VH_RANGE[0] <= vh <= VH_RANGE[1]}


def check_wing_le(cg):
    """
    Task 7 — wing_LE 補正の要否判定。

    入力: 重心 dict
    出力: dict(MAC%、補正量、補正要否)
    前提: 設計目標 CG@MAC 35%。|補正量| ≤ 1mm なら製作公差内で吸収可。
    """
    mac_pct = (cg["x_cg"] - WING_LE) / C_WING * 100.0
    correction = cg["x_cg"] - 0.35 * C_WING - WING_LE
    return {"mac_pct": mac_pct, "correction": correction,
            "need_correction": abs(correction) > 1.0}


# =============================================================================
# メイン
# =============================================================================

def _verdict(ok):
    return "✓ OK" if ok else "✗ NG"


def main():
    print("=" * 72)
    print(" ツインブーム接合方式変更による影響解析")
    print(" 2026-05-17 設計レビュー後の新方式評価")
    print("=" * 72)

    print(f"\n[新方式の前提]")
    print(f"  ブーム前端 = 主翼後縁 {WING_TE:.0f} − 埋め込み {BOOM_EMBED:.0f} = {X_BOOM_FRONT:.0f} mm")
    print(f"  ブーム後端 = 尾翼後端 = {X_BOOM_REAR:.0f} mm")
    print(f"  ブーム全長 = {BOOM_LENGTH:.0f} mm,  空中露出(片持ち)長 = {BOOM_CANTILEVER:.0f} mm")
    print(f"  主翼スパーまでのクリアランス = {X_BOOM_FRONT:.0f} − {X_WING_SPAR:.0f} = "
          f"{X_BOOM_FRONT - X_WING_SPAR:.0f} mm(干渉なし)")

    # ---- Task 1 ----
    boom = compute_boom_mass()
    print("\n" + "─" * 72)
    print(" Task 1  ブーム新質量")
    print("─" * 72)
    print(f"  m_boom = 線密度 × 全長 = {BOOM_LINEAR_DENSITY*1000:.2f} mg/mm × {BOOM_LENGTH:.0f} mm")
    print(f"         = {boom['per_boom']:.4f} g/本")
    print(f"  合計(2本)= {boom['total']:.4f} g")
    print(f"  旧値 1.8 g からの差分 Δm_boom = {boom['delta']:+.4f} g(軽量化)")

    # ---- Task 2 ----
    sleeve = compute_sleeve_mass()
    print("\n" + "─" * 72)
    print(" Task 2  EPS スリーブ質量")
    print("─" * 72)
    print(f"  V_outer = π×({SLEEVE_OD/2:.2f})²×{SLEEVE_LENGTH:.0f} = {sleeve['v_outer']:.2f} mm³")
    print(f"  V_inner = π×({SLEEVE_ID/2:.3f})²×{SLEEVE_LENGTH:.0f} = {sleeve['v_inner']:.2f} mm³")
    print(f"  V_sleeve = V_outer − V_inner = {sleeve['v_sleeve']:.2f} mm³")
    print(f"  m_sleeve = V_sleeve × {SLEEVE_DENSITY*1000:.4f} g/cm³ = {sleeve['per_sleeve']:.5f} g/個")
    print(f"  合計(2個)= {sleeve['total']:.5f} g")

    # ---- Task 3 ----
    strength = evaluate_boom_strength()
    print("\n" + "─" * 72)
    print(" Task 3  ブーム強度評価(空中露出区間、片持ち梁)")
    print("─" * 72)
    print(f"  M_boom = (L_worst/2) × ℓ_cantilever = {L_PER_BOOM:.3f} N × {BOOM_CANTILEVER/1000:.3f} m")
    print(f"         = {strength['M']*1000:.3f} mN·m")
    print(f"  W_s = π·d³/32 = {section_modulus_circular(BOOM_DIAMETER):.4f} mm³")
    print(f"  σ_boom = M / W_s = {strength['sigma']/1e6:.2f} MPa")
    print(f"  SF_boom = σ_ult / σ = 1500 / {strength['sigma']/1e6:.2f} = {strength['SF']:.0f}")
    print(f"  旧値 SF = {strength['SF_old']:.0f}(ℓ=200mm 仮定)→ 新値 SF = {strength['SF']:.0f}"
          f"(ℓ={BOOM_CANTILEVER:.0f}mm)、差分 {strength['SF']-strength['SF_old']:+.0f}")
    print(f"  判定: {_verdict(strength['SF'] >= SF_TARGET)}(目標 SF ≥ {SF_TARGET:.0f})")

    # ---- Task 4 ----
    pull = evaluate_pullout()
    rot = evaluate_rotation_mode(strength['M'])
    print("\n" + "─" * 72)
    print(" Task 4  埋め込み部 接着強度評価")
    print("─" * 72)
    print(f"  (a) 引き抜き保持力(軸方向プルアウト)")
    print(f"    A_adhesive = π×{D_HOLE}×{BOOM_EMBED:.0f} = {pull['A']:.2f} mm²")
    print(f"    τ_eff = τ×η = 10 MPa × 0.10 = {pull['tau_eff']/1e6:.1f} MPa")
    print(f"    F_pullout = τ_eff × A = {pull['F_pullout']:.2f} N")
    print(f"    F_axial(暫定)= {pull['F_axial']:.3f} N")
    print(f"    SF_pullout = F_pullout / F_axial = {pull['SF']:.0f}")
    print(f"    判定: {_verdict(pull['SF'] >= SF_TARGET)}")
    print(f"  (b) 曲げによる引き抜き(回転モード)")
    print(f"    M_embed = M_boom = {strength['M']*1000:.3f} mN·m")
    print(f"    反力作用幅 a = (2/3)×{BOOM_EMBED:.0f} = {rot['a']:.1f} mm")
    print(f"    F_couple = M_embed / a = {rot['F_couple']:.4f} N")
    print(f"    σ_eps = F_couple /(a×d_hole) = {rot['sigma_eps']/1e6:.5f} MPa")
    print(f"    SF_rotation = σ_eps_y / σ_eps = 0.4 / {rot['sigma_eps']/1e6:.5f} = {rot['SF']:.0f}")
    print(f"    判定: {_verdict(rot['SF'] >= SF_TARGET)}")

    # ---- Task 5 ----
    cg = recompute_cg(boom, sleeve)
    print("\n" + "─" * 72)
    print(" Task 5  重心位置の再計算")
    print("─" * 72)
    print(f"  新総質量 m_total = {M_TOTAL_OLD} + ({boom['delta']:+.4f}) + {sleeve['total']:.5f}")
    print(f"                   = {cg['m_total']:.4f} g")
    print(f"  Σmx 更新 = {SUM_MX_OLD} − 1.8×213 + {boom['total']:.4f}×{X_BOOM_CG:.1f}"
          f" + {sleeve['total']:.5f}×{X_SLEEVE_CG:.1f}")
    print(f"          = {cg['sum_mx']:.2f} g·mm")
    print(f"  新重心 x_cg = Σmx / m_total = {cg['x_cg']:.3f} mm")
    print(f"  旧重心 148.94 mm からのシフト = {cg['delta_x_cg']:+.3f} mm")

    # ---- Task 6 ----
    stab = check_stability(cg)
    print("\n" + "─" * 72)
    print(" Task 6  静安定性指標の確認(wing_LE = 128.0 mm 据え置き)")
    print("─" * 72)
    print(f"  中立点 x_NP/c = {stab['xnp_frac']:.4f}(Drela、AR_h={AR_H:.3f}, Vh={VH_TARGET})")
    print(f"  x_NP = {WING_LE:.0f} + {stab['xnp_frac']:.3f}×60 = {stab['x_np']:.2f} mm")
    print(f"  静安定余裕 SM = (x_NP − x_cg)/c = ({stab['x_np']:.2f} − {cg['x_cg']:.2f})/60")
    print(f"               = {stab['SM']:.1f} %  (整合済み §7 設計点 baseline {SM_BASELINE:.1f}%)")
    print(f"  判定 SM: {_verdict(stab['SM_ok'])}(目標 {SM_RANGE[0]:.0f}〜{SM_RANGE[1]:.0f} %)")
    print(f"  容積比 Vh = Sh·ℓ_h/(S·c) = {stab['Vh']:.3f}(x_cg に非依存、不変)")
    print(f"  判定 Vh: {_verdict(stab['Vh_ok'])}(目標 {VH_RANGE[0]:.2f}〜{VH_RANGE[1]:.2f})")

    # ---- Task 7 ----
    wle = check_wing_le(cg)
    print("\n" + "─" * 72)
    print(" Task 7  wing_LE 補正の要否判定")
    print("─" * 72)
    print(f"  新重心の翼弦位置 MAC% = (x_cg − wing_LE)/c = "
          f"({cg['x_cg']:.2f} − {WING_LE:.0f})/60 = {wle['mac_pct']:.2f} %")
    print(f"  wing_LE 補正量 = x_cg − 0.35·c − 128.0 = {wle['correction']:+.3f} mm")
    if wle['need_correction']:
        print(f"  判定: 補正必要(補正量 {wle['correction']:+.2f} mm、製作公差 ±1 mm 超)")
    else:
        print(f"  判定: 補正不要(|補正量| = {abs(wle['correction']):.2f} mm ≤ 1 mm、製作公差内)")

    # ---- 新旧比較表 ----
    print("\n" + "=" * 72)
    print(" 新旧比較表")
    print("=" * 72)
    rows = [
        ("ブーム質量 [g]",        f"{M_BOOM_OLD:.2f}",   f"{boom['total']:.3f}",  f"{boom['delta']:+.3f}",  "—"),
        ("スリーブ質量 [g]",      "0",                   f"{sleeve['total']:.4f}", f"+{sleeve['total']:.4f}", "—"),
        ("機体総重量 [g]",        f"{M_TOTAL_OLD:.2f}",  f"{cg['m_total']:.2f}",  f"{cg['m_total']-M_TOTAL_OLD:+.2f}", _verdict(cg['m_total'] <= 40.0)),
        ("ブーム SF",             f"{strength['SF_old']:.0f}", f"{strength['SF']:.0f}", f"{strength['SF']-strength['SF_old']:+.0f}", _verdict(strength['SF'] >= SF_TARGET)),
        ("埋込 SF(プルアウト)",  "—",                   f"{pull['SF']:.0f}",     "—",                      _verdict(pull['SF'] >= SF_TARGET)),
        ("埋込 SF(回転)",        "—",                   f"{rot['SF']:.0f}",      "—",                      _verdict(rot['SF'] >= SF_TARGET)),
        ("x_cg(設計点)[mm]",    f"{X_CG_OLD:.2f}",     f"{cg['x_cg']:.2f}",     f"{cg['delta_x_cg']:+.2f}", "—"),
        ("SM(設計点)[%]",      f"{SM_BASELINE:.1f}",  f"{stab['SM']:.1f}",     f"{stab['SM']-SM_BASELINE:+.1f}", _verdict(stab['SM_ok'])),
        ("Vh(設計点)",          f"{stab['Vh']:.3f}",   f"{stab['Vh']:.3f}",     "0",                      _verdict(stab['Vh_ok'])),
    ]
    print(f"  {'項目':<22}{'旧値':>10}{'新値':>10}{'差分':>10}  判定")
    print("  " + "-" * 64)
    for name, old, new, diff, verd in rows:
        print(f"  {name:<22}{old:>10}{new:>10}{diff:>10}  {verd}")

    print("\n  注記: 旧値の SM・Vh は 2026-05-17 v3 整合後の現行値"
          "(§7.6.2 SM=14.7%、§7.6.3 Vh=0.552)。")
    print("        ハンドオフ記載の旧値(SM 14.5%、Vh 0.567)は整合前の値であり、本表では不採用。")

    print("\n[完了]")


if __name__ == "__main__":
    main()
