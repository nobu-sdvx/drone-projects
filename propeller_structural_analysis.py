"""
propeller_structural_analysis.py
================================

§5 プロペラ設計 — 強度解析

対象: 5/18 第2回審査会の概略計算書 §5(ステップ③ 強度設計)
作成: 2026-05-16
設計根拠: Propeller_design_spec_v0.3 §5 / user-request_propeller.md §3.2

役割:
  自作 ABS プロペラ(D=36mm、3D プリント)の支配荷重(遠心力)に対する
  ブレード根応力・安全率を三点(40K/27K/22K rpm)で評価し、推定値の
  ロバスト性を RPM×ブレード質量の感度解析で確認する。併せて軸-ハブ嵌合の
  保持力と翼端マッハ数を検証する。

支配荷重:
  プロペラの破壊は推力ではなく遠心力で起こる。0.3g のブレードが自身の
  約 2,400 倍の力で根本を引っ張られる(設計書 §5.1)。

入力: propeller_config.py の凍結定数(材料 ABS、寸法、動作点)
出力(コンソール + plots/):
  - stress_distribution.png       ブレード根応力の RPM 依存
  - sf_sensitivity_heatmap.png    SF 感度ヒートマップ(m_blade × RPM)
  - hub_attachment_detail.png     軸-ハブ嵌合の詳細図

使い方:
    $ python propeller_structural_analysis.py
"""

import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["Yu Gothic", "MS Gothic", "Meiryo",
                                          "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

import propeller_config as cfg

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ブレード重心位置(設計書 §5.1、r_cg ≈ 0.75R を採用)
R_CG = 0.75 * cfg.R_PROPELLER          # = 0.0135 m
L_BLADE = cfg.R_PROPELLER - 0.002      # 揚力を受ける有効長 ≈ 16 mm


# =============================================================================
# 荷重・応力・安全率
# =============================================================================

def rpm_to_omega(rpm: float) -> float:
    """回転数 [rpm] を角速度 [rad/s] に変換する。"""
    return 2.0 * np.pi * rpm / 60.0


def compute_centrifugal_force(blade_mass: float, r_cg: float,
                              omega: float) -> float:
    """
    ブレード遠心力を計算する。

    プロペラ設計の支配荷重(推力の数百倍になる)。根本断面の引張応力評価に
    使用する。

    Args:
        blade_mass: ブレード質量 [kg]
        r_cg: ブレード重心半径 [m](通常 0.75R)
        omega: 角速度 [rad/s]

    Returns:
        遠心力 F [N]
    """
    return blade_mass * r_cg * omega ** 2


def compute_blade_root_tensile_stress(F_centrifugal: float,
                                      area_root: float) -> float:
    """根本引張応力 σ_T = F / A [Pa]。"""
    return F_centrifugal / area_root


def compute_blade_root_bending_stress(T_max: float, L_blade: float,
                                      c_root: float, t_root: float) -> float:
    """
    空気力による根本曲げ応力 σ_B = M / Z [Pa]。

    推力 T_max を有効長 L_blade に等分布荷重と仮定:
      M_root = (T/L)·L²/2 = T·L/2
      Z = b·h²/6(矩形断面、b=c_root, h=t_root)

    Args:
        T_max: ブレード 1 枚あたり最大推力 [N]
        L_blade: 揚力を受ける有効長 [m]
        c_root: 根本弦長 [m]
        t_root: 根本厚さ [m]

    Returns:
        曲げ応力 σ_B [Pa]
    """
    m_root = T_max * L_blade / 2.0
    z = c_root * t_root ** 2 / 6.0
    return m_root / z


def compute_safety_factor(sigma_combined: float,
                          sigma_y: float = cfg.SIGMA_Y_ABS) -> float:
    """安全率 SF = σ_y / σ_combined。"""
    return sigma_y / sigma_combined


def compute_hub_holding_force(adhesive_area: float,
                              tau: float = cfg.TAU_EPOXY,
                              eta: float = cfg.ETA_EFFECTIVE) -> float:
    """
    軸-ハブ嵌合の保持力 F_hold = τ · A · η [N]。

    推力反力でハブが軸から抜けるのを防ぐ。有効率 η は製作精度・接着ムラを
    織り込む保守係数(設計書 §6.3)。
    """
    return tau * adhesive_area * eta


def compute_tip_mach(omega: float, R: float,
                     c_sound: float = cfg.C_SOUND) -> float:
    """翼端マッハ数 M = ω·R / c。"""
    return omega * R / c_sound


def evaluate_stress_at_rpm(rpm: float, blade_mass: float = cfg.M_BLADE_TENTATIVE,
                           area_root: float = None) -> dict:
    """
    指定 RPM での根本応力・安全率を一括計算する。

    Returns:
        dict: 'F_cf', 'sigma_T', 'sigma_B', 'sigma_total', 'SF'
    """
    if area_root is None:
        area_root = cfg.C_ROOT * cfg.T_ROOT
    omega = rpm_to_omega(rpm)
    f_cf = compute_centrifugal_force(blade_mass, R_CG, omega)
    sigma_t = compute_blade_root_tensile_stress(f_cf, area_root)
    sigma_b = compute_blade_root_bending_stress(cfg.T_MAX_PER_MOTOR, L_BLADE,
                                                cfg.C_ROOT, cfg.T_ROOT)
    sigma_total = sigma_t + sigma_b
    return {"F_cf": f_cf, "sigma_T": sigma_t, "sigma_B": sigma_b,
            "sigma_total": sigma_total,
            "SF": compute_safety_factor(sigma_total)}


def sensitivity_analysis_rpm_mass(rpm_range: np.ndarray,
                                  mass_range: np.ndarray) -> np.ndarray:
    """
    RPM × ブレード質量の 2D 感度解析(SF 行列)を計算する。

    ブレード質量(天秤未測定)と動作 RPM(推定値)の不確実性に対し、
    安全率がどの範囲に収まるかを示す(設計書 §5.5)。

    Returns:
        SF の 2 次元配列 [mass, rpm]
    """
    area = cfg.C_ROOT * cfg.T_ROOT
    sf = np.zeros((len(mass_range), len(rpm_range)))
    for i, m in enumerate(mass_range):
        for j, rpm in enumerate(rpm_range):
            sf[i, j] = evaluate_stress_at_rpm(rpm, blade_mass=m,
                                              area_root=area)["SF"]
    return sf


# =============================================================================
# プロット
# =============================================================================

def plot_stress_distribution(save_path: str):
    """ブレード根応力・安全率の RPM 依存(計算書図 #7)。"""
    cfg.ensure_plots_dir()
    rpm = np.linspace(10000, 40000, 61)
    sig_t, sig_b, sig_tot, sf = [], [], [], []
    for r in rpm:
        res = evaluate_stress_at_rpm(r)
        sig_t.append(res["sigma_T"] / 1e6)
        sig_b.append(res["sigma_B"] / 1e6)
        sig_tot.append(res["sigma_total"] / 1e6)
        sf.append(res["SF"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))

    ax1.plot(rpm / 1000, sig_t, color="#d32f2f", lw=2, label="引張 σ_T(遠心力)")
    ax1.plot(rpm / 1000, sig_b, color="#1976d2", lw=2, label="曲げ σ_B(空気力)")
    ax1.plot(rpm / 1000, sig_tot, color="black", lw=2.4, ls="--",
             label="合成 σ_total")
    ax1.axhline(cfg.SIGMA_Y_ABS / 1e6, color="gray", ls=":",
                label=f"ABS 降伏 {cfg.SIGMA_Y_ABS / 1e6:.0f} MPa")
    ax1.set(xlabel="RPM [×1000]", ylabel="応力 [MPa]",
            title="ブレード根応力 vs RPM")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=9)

    ax2.plot(rpm / 1000, sf, color="#388e3c", lw=2.4)
    ax2.axhline(4.0, color="red", ls="--", label="目標 SF = 4")
    for name, r in cfg.RPM_EVALUATION:
        res = evaluate_stress_at_rpm(r)
        ax2.plot(r / 1000, res["SF"], "o", ms=9, color="black")
        ax2.annotate(f"{name}\nSF={res['SF']:.1f}", (r / 1000, res["SF"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax2.set(xlabel="RPM [×1000]", ylabel="安全率 SF",
            title=f"安全率 vs RPM(m_blade={cfg.M_BLADE_TENTATIVE * 1e3:.1f}g)")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=9)

    fig.suptitle("ブレード根 強度評価(支配荷重: 遠心力)", fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


def plot_sf_sensitivity_heatmap(rpm_range: np.ndarray, mass_range: np.ndarray,
                                sf_matrix: np.ndarray, save_path: str):
    """SF 感度ヒートマップ(横軸 RPM、縦軸 m_blade)。計算書図 #8。"""
    cfg.ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(8, 5))
    extent = [rpm_range[0] / 1000, rpm_range[-1] / 1000,
              mass_range[0] * 1e3, mass_range[-1] * 1e3]
    im = ax.imshow(sf_matrix, origin="lower", aspect="auto", extent=extent,
                   cmap="RdYlGn", vmin=0, vmax=15)
    cs = ax.contour(rpm_range / 1000, mass_range * 1e3, sf_matrix,
                    levels=[4, 6, 8, 10], colors="black", linewidths=1)
    ax.clabel(cs, fmt="SF=%.0f", fontsize=8)
    # SF=4 のラインを赤強調
    cs4 = ax.contour(rpm_range / 1000, mass_range * 1e3, sf_matrix,
                     levels=[4], colors="red", linewidths=2.5)
    ax.clabel(cs4, fmt="目標 SF=4", fontsize=9)
    # 三点評価をマーク
    for name, rpm in cfg.RPM_EVALUATION:
        ax.plot(rpm / 1000, cfg.M_BLADE_TENTATIVE * 1e3, "k*", ms=14)
        ax.annotate(name, (rpm / 1000, cfg.M_BLADE_TENTATIVE * 1e3),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    fig.colorbar(im, ax=ax, label="安全率 SF")
    ax.set(xlabel="RPM [×1000]", ylabel="ブレード質量 m_blade [g]",
           title="SF 感度ヒートマップ(RPM × ブレード質量)")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


def plot_hub_attachment_detail(save_path: str):
    """軸-ハブ嵌合の詳細図(計算書図 #9)。"""
    cfg.ensure_plots_dir()
    fig, ax = plt.subplots(figsize=(7, 5))

    hub_l = cfg.HUB_LENGTH * 1000
    hub_od = cfg.HUB_OUTER_DIAMETER * 1000
    hole_d = cfg.HUB_HOLE_DIAMETER * 1000
    shaft_d = cfg.SHAFT_DIAMETER * 1000

    # ハブ本体(断面)
    ax.add_patch(plt.Rectangle((0, -hub_od / 2), hub_l, hub_od,
                               facecolor="#bcd4e6", edgecolor="black", lw=1.5,
                               label="ABS ハブ"))
    # 接着層(誇張表示)
    ax.add_patch(plt.Rectangle((0, -hole_d / 2 - 0.15), hub_l, 0.15,
                               facecolor="#f57c00", edgecolor="none"))
    ax.add_patch(plt.Rectangle((0, hole_d / 2), hub_l, 0.15,
                               facecolor="#f57c00", edgecolor="none",
                               label="エポキシ接着層"))
    # モーター軸
    ax.add_patch(plt.Rectangle((-3, -shaft_d / 2), hub_l + 3, shaft_d,
                               facecolor="gray", edgecolor="black",
                               label="モーター軸 φ1.0mm"))
    # 寸法注記
    ax.annotate("", xy=(0, hub_od / 2 + 1), xytext=(hub_l, hub_od / 2 + 1),
                arrowprops=dict(arrowstyle="<->", color="blue"))
    ax.text(hub_l / 2, hub_od / 2 + 1.4, f"ハブ長 {hub_l:.0f} mm",
            ha="center", color="blue", fontsize=9)
    ax.text(hub_l + 0.5, 0, f"穴 φ{hole_d:.2f}mm\n(締めしろ -0.05mm)",
            fontsize=8, va="center")

    ax.set(xlim=(-5, hub_l + 8), ylim=(-hub_od, hub_od + 3),
           xlabel="軸方向 [mm]", ylabel="径方向 [mm]",
           title="軸-ハブ嵌合詳細(圧入 + エポキシ接着)")
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → 保存: {save_path}")


# =============================================================================
# メイン
# =============================================================================

def main():
    print("=" * 72)
    print(" §5 プロペラ設計  強度解析  (作成 2026-05-16)")
    print("=" * 72)

    area_root = cfg.C_ROOT * cfg.T_ROOT
    print(f"\n[強度設計の前提]")
    print(f"  材料            = ABS(3Dプリント、σ_y={cfg.SIGMA_Y_ABS / 1e6:.0f} MPa)")
    print(f"  根本断面        = {cfg.C_ROOT * 1e3:.1f} × {cfg.T_ROOT * 1e3:.1f} mm"
          f" = {area_root * 1e6:.2f} mm²")
    print(f"  ブレード質量    = {cfg.M_BLADE_TENTATIVE * 1e3:.1f} g(暫定、天秤実測待ち)")
    print(f"  重心半径 r_cg   = {R_CG * 1e3:.1f} mm(0.75R)")

    # ---------- 三点評価 ----------
    print("\n" + "─" * 72)
    print(" 三点評価(40K最悪 / 27K最大動作 / 22K巡航)")
    print("─" * 72)
    print(f"\n  {'動作点':<14}{'F_cf [N]':>10}{'σ_T [MPa]':>11}"
          f"{'σ_B [MPa]':>11}{'σ_total':>10}{'SF':>8}")
    results = {}
    for name, rpm in cfg.RPM_EVALUATION:
        r = evaluate_stress_at_rpm(rpm)
        results[name] = r
        flag = "✓" if r["SF"] >= 4 else "△ 目標4未満"
        print(f"  {name:<14}{r['F_cf']:>10.1f}{r['sigma_T'] / 1e6:>11.2f}"
              f"{r['sigma_B'] / 1e6:>11.2f}{r['sigma_total'] / 1e6:>10.2f}"
              f"{r['SF']:>8.1f}  {flag}")
    print(f"\n  例え話: 40K rpm の遠心力 {results['40K最悪']['F_cf']:.0f} N は"
          f" 約 {results['40K最悪']['F_cf'] / cfg.G:.1f} kg重。")
    print(f"          {cfg.M_BLADE_TENTATIVE * 1e3:.1f}g のブレードが自身の約"
          f" {results['40K最悪']['F_cf'] / (cfg.M_BLADE_TENTATIVE * cfg.G):.0f}"
          f" 倍の力で根本を引かれる。")

    # ---------- 翼端マッハ数 ----------
    print("\n" + "─" * 72)
    print(" 翼端マッハ数")
    print("─" * 72)
    for name, rpm in cfg.RPM_EVALUATION:
        m_tip = compute_tip_mach(rpm_to_omega(rpm), cfg.R_PROPELLER)
        v_tip = rpm_to_omega(rpm) * cfg.R_PROPELLER
        print(f"  {name:<14} V_tip = {v_tip:5.1f} m/s,  M_tip = {m_tip:.3f}"
              f"  {'✓ 圧縮性影響なし' if m_tip < 0.3 else '△'}")

    # ---------- 軸-ハブ嵌合 ----------
    print("\n" + "─" * 72)
    print(" 軸-ハブ嵌合の保持力")
    print("─" * 72)
    a_adh = np.pi * cfg.SHAFT_DIAMETER * cfg.HUB_LENGTH
    f_hold = compute_hub_holding_force(a_adh)
    sf_hub = f_hold / cfg.T_MAX_PER_MOTOR
    print(f"  接着面積 A = π·d·L = π × {cfg.SHAFT_DIAMETER * 1e3:.1f}mm ×"
          f" {cfg.HUB_LENGTH * 1e3:.0f}mm = {a_adh * 1e6:.2f} mm²")
    print(f"  保持力 F_hold = τ·A·η = {cfg.TAU_EPOXY / 1e6:.0f}MPa ×"
          f" {a_adh * 1e6:.2f}mm² × {cfg.ETA_EFFECTIVE:.2f}(有効率)"
          f" = {f_hold:.1f} N")
    print(f"  推力反力 = {cfg.T_MAX_PER_MOTOR * 1e3:.0f} mN")
    print(f"  嵌合 SF = F_hold / T_max = {sf_hub:.0f}"
          f"  {'✓ 大幅余裕' if sf_hub >= 100 else '△'}")

    # ---------- 感度解析 ----------
    print("\n" + "─" * 72)
    print(" 感度解析(ブレード質量 × RPM)")
    print("─" * 72)
    rpm_range = np.linspace(15000, 40000, 26)
    mass_range = np.linspace(0.2e-3, 0.5e-3, 7)
    sf_matrix = sensitivity_analysis_rpm_mass(rpm_range, mass_range)
    print(f"\n  {'m_blade[g]':<12}", end="")
    for name, rpm in cfg.RPM_EVALUATION:
        print(f"{'SF@' + name:>14}", end="")
    print()
    for m in mass_range:
        print(f"  {m * 1e3:<12.2f}", end="")
        for name, rpm in cfg.RPM_EVALUATION:
            sf = evaluate_stress_at_rpm(rpm, blade_mass=m,
                                        area_root=area_root)["SF"]
            print(f"{sf:>14.1f}", end="")
        print()
    print(f"\n  → m_blade=0.4g(最悪想定)でも 27K 動作で SF="
          f"{evaluate_stress_at_rpm(cfg.RPM_MAX_OP, blade_mass=0.4e-3)['SF']:.1f}"
          f" を維持。")

    # ---------- プロット ----------
    print("\n" + "─" * 72)
    print(" 図表生成(plots/)")
    print("─" * 72)
    plot_stress_distribution(cfg.plot_path("stress_distribution.png"))
    plot_sf_sensitivity_heatmap(rpm_range, mass_range, sf_matrix,
                                cfg.plot_path("sf_sensitivity_heatmap.png"))
    plot_hub_attachment_detail(cfg.plot_path("hub_attachment_detail.png"))

    print("\n" + "=" * 72)
    print(" 確定値サマリー(計算書 §5 提供物)")
    print("=" * 72)
    print(f"  ブレード根遠心力(40K)     = {results['40K最悪']['F_cf']:.1f} N")
    print(f"  合成応力(40K最悪)         = "
          f"{results['40K最悪']['sigma_total'] / 1e6:.2f} MPa")
    print(f"  安全率 SF  40K / 27K / 22K  = "
          f"{results['40K最悪']['SF']:.1f} / {results['27K最大動作']['SF']:.1f}"
          f" / {results['22K巡航']['SF']:.1f}")
    print(f"  翼端マッハ数(40K)         = "
          f"{compute_tip_mach(rpm_to_omega(cfg.RPM_MAX), cfg.R_PROPELLER):.3f}")
    print(f"  軸-ハブ嵌合 SF              = {sf_hub:.0f}")
    print("\n[完了]")

    return {"results": results, "sf_hub": sf_hub}


if __name__ == "__main__":
    main()
