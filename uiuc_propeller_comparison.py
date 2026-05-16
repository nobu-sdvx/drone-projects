"""
uiuc_propeller_comparison.py
============================

§5 プロペラ設計 — UIUC 実測データとの交差検証

対象: 5/18 第2回審査会の概略計算書 §5
作成: 2026-05-16
設計根拠: Propeller_design_spec_v0.3 §4.7 / user-request_propeller.md §3.3, §4

役割:
  自作プロペラの形状・推力予測を、UIUC Propeller Database Volume 2 の
  実測プロペラ(主参考: Crazyflie cfnq_45p1 = 45mm)と比較する。
  自作 BEM/QPROP は低 Re 域で ±20% の不確実性があるため(設計書 §9.3)、
  実測データとの照合で「予測トレンドが実測と整合」することを示す。

入力(uiuc_propeller_data/volume-2/data/):
  cfnq_45p1_geom.txt            … Crazyflie 45mm ブレード形状
  cfnq_45p1_static_1322rd.txt   … Crazyflie 45mm 静推力試験データ

出力(コンソール + plots/):
  - crazyflie_comparison_geom.png    形状比較(自作 36mm vs Crazyflie 45mm)
  - crazyflie_comparison_thrust.png  推力比較(自作 BEM/QPROP vs Crazyflie 実測)

使い方:
    $ python uiuc_propeller_comparison.py
"""

import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["Yu Gothic", "MS Gothic", "Meiryo",
                                          "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

import propeller_config as cfg
from propeller_design_analysis import (load_blade_geometry_from_crazyflie,
                                       run_qprop, screen_airfoil_candidates,
                                       bem_analysis)

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

D_CRAZYFLIE = 0.045   # Crazyflie 45mm プロペラ直径 [m]


# =============================================================================
# UIUC データ読み込み
# =============================================================================

def load_uiuc_geometry(filepath: str) -> dict:
    """
    UIUC プロペラ geometry ファイル(.txt)を読み込む。

    フォーマット:
        r/R   c/R     beta
        0.15  0.2453  29.162
        ...

    Returns:
        dict: {'r_R', 'c_R', 'beta_deg'}(いずれも np.ndarray)
    """
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append([float(p) for p in parts[:3]])
            except ValueError:
                continue   # ヘッダ行をスキップ
    data = np.array(rows)
    return {"r_R": data[:, 0], "c_R": data[:, 1], "beta_deg": data[:, 2]}


def load_uiuc_static_data(filepath: str) -> dict:
    """
    UIUC 静推力試験データ(.txt)を読み込む。

    フォーマット:
        RPM        CT        CP
         4026.667  0.139894  0.131208
        ...
    静推力 = advance ratio J=0(機体速度ゼロ)での測定。

    Returns:
        dict: {'rpm', 'CT', 'CP'}(いずれも np.ndarray)
    """
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append([float(p) for p in parts[:3]])
            except ValueError:
                continue
    data = np.array(rows)
    return {"rpm": data[:, 0], "CT": data[:, 1], "CP": data[:, 2]}


def ct_to_thrust(ct, rpm, D: float, rho: float = cfg.RHO):
    """
    推力係数 CT から推力 T を計算する。

    T = CT · ρ · n² · D⁴,  n = RPM / 60 [rev/s]

    Args:
        ct: 推力係数(スカラまたは配列)
        rpm: 回転数 [rpm]
        D: プロペラ直径 [m]

    Returns:
        推力 T [N]
    """
    n = np.asarray(rpm) / 60.0
    return ct * rho * n ** 2 * D ** 4


# =============================================================================
# プロット
# =============================================================================

def compare_geometry_with_uiuc(self_geom, uiuc_geom: dict, save_path: str):
    """自作 36mm と Crazyflie 45mm のブレード形状比較(計算書図 #10)。"""
    cfg.ensure_plots_dir()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))

    # 弦長(原寸 mm)
    ax1.plot(self_geom.radial_stations, self_geom.chord_distribution * 1e3,
             "o-", color="#d32f2f", lw=2, label="自作 D=36mm")
    ax1.plot(uiuc_geom["r_R"], uiuc_geom["c_R"] * D_CRAZYFLIE / 2 * 1e3,
             "s--", color="#1976d2", lw=1.6, label="Crazyflie D=45mm")
    ax1.set(xlabel="r/R", ylabel="chord c [mm]", title="弦長分布の比較")
    ax1.grid(alpha=0.3)
    ax1.legend()

    # ねじり角(自作は Crazyflie をそのまま流用しているので一致)
    ax2.plot(self_geom.radial_stations,
             np.degrees(self_geom.twist_distribution),
             "o-", color="#d32f2f", lw=2, label="自作 D=36mm")
    ax2.plot(uiuc_geom["r_R"], uiuc_geom["beta_deg"],
             "s--", color="#1976d2", lw=1.6, label="Crazyflie cfnq_45p1")
    ax2.set(xlabel="r/R", ylabel="twist β [deg]", title="ねじり角分布の比較")
    ax2.grid(alpha=0.3)
    ax2.legend()

    fig.suptitle("ブレード形状: 自作 36mm vs Crazyflie 45mm(UIUC Vol.2)",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


def compare_thrust_curves_with_uiuc(cf_static: dict, self_static: dict,
                                    scaled_pred: dict, save_path: str):
    """
    推力比較プロット(計算書図 #11)。

    - Crazyflie 45mm: UIUC 実測静推力
    - 自作 36mm: QPROP 静推力(V=0)
    - 参考: Crazyflie 実測を D・RPM スケーリングした単純予測(設計書 §4.7)
    """
    cfg.ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(7.8, 5))

    t_cf = ct_to_thrust(cf_static["CT"], cf_static["rpm"], D_CRAZYFLIE)
    ax.plot(cf_static["rpm"] / 1000, t_cf * 1e3, "s-", color="#1976d2",
            lw=1.8, label="Crazyflie 45mm(UIUC 実測静推力)")
    ax.plot(self_static["rpm"] / 1000, self_static["thrust"] * 1e3, "o-",
            color="#d32f2f", lw=2, label="自作 36mm(QPROP 静推力)")
    ax.plot(scaled_pred["rpm"] / 1000, scaled_pred["thrust"] * 1e3, "^:",
            color="#f57c00", lw=1.6,
            label="自作 36mm(Crazyflie 実測の単純スケーリング)")
    ax.axhline(cfg.T_REQUIRED_PER_MOTOR * 1e3, color="green", ls="--",
               label=f"巡航必要推力 {cfg.T_REQUIRED_PER_MOTOR * 1e3:.1f} mN")
    ax.set(xlabel="RPM [×1000]", ylabel="静推力 T [mN]",
           title="静推力比較: 自作 vs Crazyflie 実測")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


# =============================================================================
# メイン
# =============================================================================

def main():
    print("=" * 72)
    print(" §5 プロペラ設計  UIUC 実測データとの交差検証  (作成 2026-05-16)")
    print("=" * 72)

    # ---------- UIUC データ読み込み ----------
    cf_geom = load_uiuc_geometry(cfg.CRAZYFLIE_GEOM_FILE)
    cf_static = load_uiuc_static_data(cfg.CRAZYFLIE_STATIC_FILE)
    print(f"\n[UIUC データ読み込み(Crazyflie cfnq_45p1)]")
    print(f"  geometry   : {len(cf_geom['r_R'])} ステーション "
          f"(r/R = {cf_geom['r_R'][0]:.2f}〜{cf_geom['r_R'][-1]:.2f})")
    print(f"  静推力試験 : {len(cf_static['rpm'])} 点 "
          f"(RPM = {cf_static['rpm'][0]:.0f}〜{cf_static['rpm'][-1]:.0f})")

    # ---------- 自作プロペラ形状(翼型は設計スクリプトで選定) ----------
    _, best_airfoil = screen_airfoil_candidates()
    self_geom = load_blade_geometry_from_crazyflie(D=cfg.D_PROPELLER,
                                                   airfoil_name=best_airfoil)

    # ---------- 形状比較 ----------
    print("\n" + "─" * 72)
    print(" ブレード形状の比較(75%R 代表点)")
    print("─" * 72)
    i75_self = int(np.argmin(np.abs(self_geom.radial_stations - 0.75)))
    i75_cf = int(np.argmin(np.abs(cf_geom["r_R"] - 0.75)))
    c75_self = self_geom.chord_distribution[i75_self] * 1e3
    c75_cf = cf_geom["c_R"][i75_cf] * D_CRAZYFLIE / 2 * 1e3
    b75_self = np.degrees(self_geom.twist_distribution[i75_self])
    b75_cf = cf_geom["beta_deg"][i75_cf]
    print(f"  {'項目':<16}{'自作 36mm':>14}{'Crazyflie 45mm':>16}{'差':>10}")
    print(f"  {'弦長 c75 [mm]':<16}{c75_self:>14.2f}{c75_cf:>16.2f}"
          f"{(c75_self - c75_cf) / c75_cf * 100:>9.1f}%")
    print(f"  {'ねじり角 β75 [°]':<15}{b75_self:>14.2f}{b75_cf:>16.2f}"
          f"{b75_self - b75_cf:>9.2f}°")
    print(f"  → c/R, β は無次元形状をそのまま流用。弦長は D 比 36/45 で縮小。")

    # ---------- Crazyflie 実測の代表点 ----------
    print("\n" + "─" * 72)
    print(" Crazyflie 45mm 実測静推力(最高 RPM 点)")
    print("─" * 72)
    rpm_cf = cf_static["rpm"][-1]
    ct_cf = cf_static["CT"][-1]
    t_cf = ct_to_thrust(ct_cf, rpm_cf, D_CRAZYFLIE)
    print(f"  RPM = {rpm_cf:.0f},  CT = {ct_cf:.4f}")
    print(f"  推力 T = CT·ρ·n²·D⁴ = {t_cf * 1e3:.2f} mN "
          f"({t_cf / cfg.G * 1e3:.2f} g重)")

    # ---------- 設計書 §4.7 単純スケーリング検算 ----------
    print("\n" + "─" * 72)
    print(" 推力スケーリング検算(設計書 §4.7)")
    print("─" * 72)
    # T ∝ ρ·n²·D⁴(同一形状) → Crazyflie 実測を自作動作点に換算
    n_self = cfg.RPM_CRUISE / 60.0
    n_cf = rpm_cf / 60.0
    ratio = (n_self / n_cf) ** 2 * (cfg.D_PROPELLER / D_CRAZYFLIE) ** 4
    t_scaled = t_cf * ratio
    print(f"  T_自作/T_Crazyflie = (n_自作/n_Crazy)²·(D_自作/D_Crazy)⁴")
    print(f"                    = ({n_self:.0f}/{n_cf:.0f})² ×"
          f" ({cfg.D_PROPELLER * 1e3:.0f}/{D_CRAZYFLIE * 1e3:.0f})⁴ = {ratio:.3f}")
    print(f"  → 単純スケーリング推力 = {t_scaled * 1e3:.2f} mN "
          f"({t_scaled / cfg.G * 1e3:.2f} g重)")
    print(f"  ※ これは静推力(J=0)ベース。固定ピッチ・固定RPMでは前進飛行"
          f"(前進比 J↑)で推力は減少する。")

    # ---------- 自作プロペラの QPROP 静推力スイープ ----------
    print("\n" + "─" * 72)
    print(" 自作 36mm の静推力(QPROP、V=0)")
    print("─" * 72)
    refined_prop = os.path.join(cfg.QPROP_DATA_DIR,
                                "powerup_propeller_refined.prop")
    if not os.path.exists(refined_prop):
        refined_prop = os.path.join(cfg.QPROP_DATA_DIR,
                                    "powerup_propeller.prop")
        print(f"  [注意] refined.prop 未生成 → 暫定 .prop を使用"
              f"(先に propeller_design_analysis.py を実行推奨)")
    motor_file = os.path.join(cfg.QPROP_DATA_DIR, "powerup_motor.mot")
    rpm_self = np.array([10000, 15000, 20000, 22100, 26600, 30000,
                         35000, 40000], dtype=float)
    t_self = np.array([run_qprop(refined_prop, motor_file, 0.0, r)["thrust"]
                       for r in rpm_self])
    self_static = {"rpm": rpm_self, "thrust": t_self}
    for r, t in zip(rpm_self, t_self):
        print(f"  {r:7.0f} rpm : T = {t * 1e3:7.2f} mN")

    # スケーリング予測曲線(Crazyflie 実測を D 比で換算)
    scaled_pred = {"rpm": cf_static["rpm"],
                   "thrust": ct_to_thrust(cf_static["CT"], cf_static["rpm"],
                                          cfg.D_PROPELLER)}

    # ---------- 同一 RPM での照合 ----------
    print("\n" + "─" * 72)
    print(" 同一 RPM(20,000)での自作 vs Crazyflie 照合")
    print("─" * 72)
    t_self_20k = run_qprop(refined_prop, motor_file, 0.0, 20000)["thrust"]
    ct_cf_20k = float(np.interp(20000, cf_static["rpm"], cf_static["CT"]))
    t_cf_20k_at36 = ct_to_thrust(ct_cf_20k, 20000, cfg.D_PROPELLER)
    print(f"  自作 36mm QPROP            : T = {t_self_20k * 1e3:.2f} mN")
    print(f"  Crazyflie CT を D=36mm 換算 : T = {t_cf_20k_at36 * 1e3:.2f} mN")
    dev = abs(t_self_20k - t_cf_20k_at36) / t_cf_20k_at36 * 100
    print(f"  乖離 = {dev:.1f} %  "
          f"({'✓ 実測トレンドと整合' if dev <= 35 else '△ 要考察'})")

    # ---------- プロット ----------
    print("\n" + "─" * 72)
    print(" 図表生成(plots/)")
    print("─" * 72)
    compare_geometry_with_uiuc(self_geom, cf_geom,
                               cfg.plot_path("crazyflie_comparison_geom.png"))
    compare_thrust_curves_with_uiuc(cf_static, self_static, scaled_pred,
                                    cfg.plot_path("crazyflie_comparison_thrust.png"))

    print("\n" + "=" * 72)
    print(" 確定値サマリー(計算書 §5 提供物)")
    print("=" * 72)
    print(f"  Crazyflie 45mm @ {rpm_cf:.0f}rpm 実測 : T = {t_cf * 1e3:.1f} mN")
    print(f"  自作 36mm 単純スケーリング推力       : T = {t_scaled * 1e3:.1f} mN")
    print(f"  自作 36mm QPROP 静推力(22.1K rpm)  : "
          f"T = {self_static['thrust'][3] * 1e3:.1f} mN")
    print(f"  75%R ねじり角差(自作 vs Crazyflie)  : {b75_self - b75_cf:+.2f}°")
    print("\n[完了]")

    return {"cf_geom": cf_geom, "cf_static": cf_static,
            "self_static": self_static, "t_scaled": t_scaled}


if __name__ == "__main__":
    main()
