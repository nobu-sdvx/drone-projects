"""
weight_cg_analysis.py
======================
§7 重量配分・重心解析

応用機械設計製図Ⅰ ドローン班 (2026)
作成日: 2026-05-06
作成者: Nobu (Claude支援)

目的:
  POWERUP 4.0 (実測17g) と各部品の質量・配置から、機体総重量と
  重心位置 (x_cg) を計算し、目標CG位置(MAC 30%、35%)を満たす
  主翼前縁 (wing_LE) 位置を逆算する。

入力:
  - 部品の質量・x座標(機首=0、後方が正)
  - CG目標位置(MAC比)
  - 機体総重量シナリオ(楽観/主/悲観)

出力:
  - 主翼LE位置 [mm]
  - 機体総重量 [g]
  - 静安定余裕 SM [% MAC]
  - 尾翼モーメントアーム ℓh [mm]
  - 容積比 Vh, Vv
  - 機体側面図 PNG
  - CGスイーププロット PNG

依存:
  pip install numpy matplotlib japanize-matplotlib
"""

import os
from dataclasses import dataclass
from typing import List, Dict
import numpy as np
import matplotlib.pyplot as plt

# 日本語フォント(なければ英語にフォールバック)
try:
    import japanize_matplotlib  # noqa
    JP_OK = True
except ImportError:
    JP_OK = False
    print("[WARN] japanize-matplotlib 未インストール。ラベルは英語にフォールバック。")


# =============================================================================
# Section 1: 機体定数(根拠付き)
# =============================================================================
# ★v3 更新(2026-05-17):ツインブーム + U字型尾翼 + ロッドフェアリング構成。
#   §3 尾翼U字版v3・§5 プロペラ・§7 重量重心U字版 と数値を統一。

# --- 機体寸法(§2 主翼 概略計算書)---
MAC_MM = 60.0           # 平均空力翼弦 [mm](矩形翼の翼弦)
WING_AREA_CM2 = 216.0   # 翼面積 [cm²](360×60mm)
SPAN_MM = 360.0         # 翼幅 [mm]
AR = 6.0                # アスペクト比

# --- 尾翼寸法(§3 尾翼U字版v3 で確定)---
SH_CM2 = 53.28          # 水平尾翼面積 [cm²](148×36mm)
SV_CM2 = 24.0           # 垂直尾翼面積 [cm²](40×30mm × U字2枚)
AR_H = 148.0 / 36.0     # 水平尾翼アスペクト比 = 4.111(148×36mm)

# 機体最後端 305 mm。水平尾翼AC = 後端 − 0.75·c_h。
# v3 では尾翼AC位置 277.25 mm を機体配置のアンカーとして固定する
# (§3 v3 改訂:c_h を 37→36 mm に微調整しても尾翼AC位置は不変)。
X_TAIL_END_MM = 305.0   # 機体最後端 = 尾翼質量の代表位置 [mm]
X_TAIL_AC_MM = 277.25   # 水平尾翼AC位置 [mm](§3 v3 固定アンカー)

# --- 中立点位置(Drela MIT 16.01 "Lab 8 Notes" 簡易式)---
# x_np/MAC = 1/4 + (1 - 4/(AR+2)) × Vh × (1+2/AR)/(1+2/AR_h)
VH_TARGET = 0.55        # 水平尾翼容積比 設計目標(§3)


def neutral_point_frac(ar=AR, ar_h=AR_H, vh=VH_TARGET):
    """中立点位置 x_np/MAC を Drela 近似式で計算する。"""
    return 0.25 + (1.0 - 4.0 / (ar + 2.0)) * vh * (1.0 + 2.0 / ar) / (1.0 + 2.0 / ar_h)


X_NP_MAC_FRAC = neutral_point_frac()   # v3: AR_h=4.111 → 0.497

# --- POWERUP 4.0(★Nobu実測値 2026-05-06★)---
POWERUP_MASS_G = 17.0       # 実測(プロペラ抜き、公式19gはストック品込み)
POWERUP_CG_X_MM = 111.0     # 実測(機首=POWERUP前端から)
POWERUP_LENGTH_MM = 220.0   # 公式仕様

# --- CG目標 ---
CG_TARGETS = [0.30, 0.35]   # 30%、35% MAC

# --- 機体総重量シナリオ(感度解析用) ---
# 主シナリオ(36.5g)を基準に、±2g程度の不確実性を見込む
WEIGHT_SCENARIOS = {
    "optimistic": {"misc_extra": -1.5, "label": "楽観 (W_min)"},
    "main":       {"misc_extra":  0.0, "label": "主シナリオ"},
    "pessimistic":{"misc_extra": +2.0, "label": "悲観 (W_max)"},
}


# =============================================================================
# Section 2: 部品データ構造
# =============================================================================

@dataclass
class Component:
    """機体部品の質量と位置"""
    name: str
    mass_g: float
    x_mm: float          # 機首=0、後方が正
    notes: str = ""


def get_components(wing_le_mm: float, scenario: str = "main") -> List[Component]:
    """
    配置案に基づく部品リストを返す。

    Args:
        wing_le_mm: 主翼前縁の機首からの距離 [mm]
        scenario: 重量シナリオ ("optimistic"/"main"/"pessimistic")

    Returns:
        Component のリスト

    根拠:
        - 質量: requirements.md §3.4 改訂版 + POWERUP実測
        - x座標: §7 配置案(プッシャー、POWERUP前向き)
    """
    misc_extra = WEIGHT_SCENARIOS[scenario]["misc_extra"]
    misc_mass = max(0.5, 3.0 + misc_extra)  # 接着剤等は0.5g下限

    # v3 部品リスト(§7.2 ツインブーム+U字+ロッドフェアリング構成、全11部品)
    return [
        Component("POWERUP 4.0",            POWERUP_MASS_G, POWERUP_CG_X_MM,
                  "実測値、プロペラ抜き"),
        Component("胴体取付フェアリング",     3.0, 110.0,
                  "POWERUP固定用ジグ・サドルバンド"),
        Component("ロッドフェアリング本体",   0.70, 132.7,
                  "EPSかまぼこ型、露出ロッド被覆"),
        Component("主翼サドル",              0.05, 158.3,
                  "フェアリング上面の主翼着座凸部"),
        Component("自作プロペラ+マウント",   5.0, 230.0,
                  "POWERUPモーター(x≈220)直後、2基分"),
        Component("主翼 AG36",               5.0, wing_le_mm + MAC_MM/2,
                  "矩形翼、幾何中心 = LE + MAC/2"),
        Component("ツインブーム",            1.432, 231.5,
                  "φ2mm CFRP × 2本、主翼後縁から30mm埋め込み"),
        Component("EPS スリーブ ×2",         0.012, 195.5,
                  "ブーム埋め込み部の根元補強(OD φ5 × L 15mm)"),
        Component("水平尾翼",                0.6, X_TAIL_END_MM,
                  "EPS+φ1mm CFRPスパー、機体最後端"),
        Component("垂直尾翼 ×2(U字)",       0.32, X_TAIL_END_MM,
                  "0.16g/枚 × 2、ブーム後端から上向き"),
        Component("接着剤・補強・配線",       misc_mass, 170.0,
                  f"機体中央付近に分布(scenario={scenario})"),
    ]


# =============================================================================
# Section 3: 計算関数
# =============================================================================

def compute_cg(components: List[Component]) -> float:
    """重心位置 x_cg [mm] を計算: x_cg = Σ(m_i × x_i) / Σ(m_i)"""
    total_mass = sum(c.mass_g for c in components)
    total_moment = sum(c.mass_g * c.x_mm for c in components)
    return total_moment / total_mass


def compute_total_mass(components: List[Component]) -> float:
    return sum(c.mass_g for c in components)


def find_wing_le_for_cg_target(target_mac_frac: float,
                                scenario: str = "main") -> float:
    """
    目標CG位置(MAC%)を満たす主翼LE位置を解析的に解く。

    式の導出:
        x_cg_target = wing_LE + target_mac_frac × MAC
        x_cg_actual = [m_wing × (wing_LE + MAC/2) + Σ(他部品 m_i × x_i)] / W

        この2つを等しくして wing_LE について解く。

    Args:
        target_mac_frac: 目標CG位置(MAC比)、例: 0.30 = 30% MAC

    Returns:
        wing_LE [mm]
    """
    # 主翼以外を取得して、固定モーメントと固定質量を計算
    components_with_dummy = get_components(wing_le_mm=0.0, scenario=scenario)
    fixed = [c for c in components_with_dummy if c.name != "主翼 AG36"]
    fixed_mass = sum(c.mass_g for c in fixed)
    fixed_moment = sum(c.mass_g * c.x_mm for c in fixed)

    wing_mass = next(c.mass_g for c in components_with_dummy if c.name == "主翼 AG36")
    total_mass = fixed_mass + wing_mass

    # 連立方程式を解く
    # total_mass × (wing_LE + t·MAC) = wing_mass × (wing_LE + MAC/2) + fixed_moment
    # → fixed_mass × wing_LE = wing_mass × MAC/2 + fixed_moment - total_mass × t × MAC
    wing_le = (wing_mass * MAC_MM / 2 + fixed_moment
               - total_mass * target_mac_frac * MAC_MM) / fixed_mass
    return wing_le


def compute_static_margin(wing_le_mm: float, x_cg_mm: float) -> Dict:
    """
    静安定余裕 SM [% MAC] を計算

    SM = (x_np - x_cg) / MAC

    SM > 0: 縦に静安定(自分で水平に戻る性質あり)
    SM < 0: 不安定(操縦不能)
    SM 5-15%: 一般的な目標範囲(動きやすく安定)
    SM > 20%: 過安定(動きが鈍い)
    """
    x_np_mm = wing_le_mm + X_NP_MAC_FRAC * MAC_MM
    sm_frac = (x_np_mm - x_cg_mm) / MAC_MM
    return {
        "x_np_mm": x_np_mm,
        "SM_pct": sm_frac * 100,
    }


def compute_tail_volume(wing_le_mm: float) -> Dict:
    """
    尾翼容積比 Vh, Vv を計算
        Vh = Sh × ℓh / (S × MAC)
        Vv = Sv × ℓv / (S × b)

    Drela 推奨範囲: Vh = 0.30 ~ 0.60、Vv = 0.02 ~ 0.05
    """
    wing_ac_mm = wing_le_mm + 0.25 * MAC_MM
    lh_mm = X_TAIL_AC_MM - wing_ac_mm

    # 単位を [cm] に揃える
    lh_cm = lh_mm / 10.0
    mac_cm = MAC_MM / 10.0
    span_cm = SPAN_MM / 10.0

    vh = (SH_CM2 * lh_cm) / (WING_AREA_CM2 * mac_cm)
    vv = (SV_CM2 * lh_cm) / (WING_AREA_CM2 * span_cm)

    return {
        "wing_ac_mm": wing_ac_mm,
        "lh_mm": lh_mm,
        "Vh": vh,
        "Vv": vv,
    }


# =============================================================================
# Section 4: 可視化
# =============================================================================

def plot_side_view(wing_le_mm: float, scenario: str, save_path: str,
                   target_label: str):
    """機体側面図 + 部品配置 + CG位置"""
    components = get_components(wing_le_mm, scenario)
    x_cg = compute_cg(components)
    total_mass = compute_total_mass(components)
    tail = compute_tail_volume(wing_le_mm)
    sm = compute_static_margin(wing_le_mm, x_cg)

    fig, ax = plt.subplots(figsize=(13, 5))

    # 胴体ロッドを表現(POWERUP+延長ブーム)
    ax.plot([0, X_TAIL_END_MM], [0, 0], '-', color='gray', lw=4, alpha=0.4)
    ax.plot([0, POWERUP_LENGTH_MM], [0, 0], '-', color='black', lw=3, alpha=0.7)

    # 主翼を矩形で表現
    wing_rect_x = [wing_le_mm, wing_le_mm + MAC_MM,
                   wing_le_mm + MAC_MM, wing_le_mm, wing_le_mm]
    wing_rect_y = [-30, -30, 30, 30, -30]
    ax.fill(wing_rect_x, wing_rect_y, alpha=0.25, color='steelblue',
            label='主翼 AG36' if JP_OK else 'Wing AG36')

    # 尾翼を矩形で
    ax.fill([X_TAIL_END_MM-15, X_TAIL_END_MM+5, X_TAIL_END_MM+5, X_TAIL_END_MM-15, X_TAIL_END_MM-15],
            [-15, -15, 15, 15, -15], alpha=0.25, color='darkorange',
            label='尾翼' if JP_OK else 'Tail')

    # 各部品を縦バーでマーク
    colors = plt.cm.tab10(np.linspace(0, 1, len(components)))
    y_offset = 50
    for c, color in zip(components, colors):
        bar_h = c.mass_g * 5  # 質量を高さで表現
        ax.bar(c.x_mm, bar_h, width=4, bottom=y_offset, color=color,
               alpha=0.8, edgecolor='black', lw=0.5)
        ax.text(c.x_mm, y_offset + bar_h + 5, f"{c.mass_g}g",
                ha='center', fontsize=8)
        ax.text(c.x_mm, y_offset - 8, c.name if JP_OK else c.name.encode('ascii', 'replace').decode(),
                ha='center', fontsize=7, rotation=45)

    # CG位置をマーク
    ax.axvline(x_cg, color='red', linestyle='--', lw=2,
               label=f'CG = {x_cg:.1f}mm' if JP_OK else f'CG = {x_cg:.1f}mm')
    ax.plot(x_cg, 0, 'rv', markersize=15)

    # 中立点もマーク
    ax.axvline(sm['x_np_mm'], color='green', linestyle=':', lw=1.5,
               label=f'NP = {sm["x_np_mm"]:.1f}mm' if JP_OK else f'NP = {sm["x_np_mm"]:.1f}mm')

    # 設定
    ax.set_xlim(-20, X_TAIL_END_MM + 30)
    ax.set_ylim(-50, 200)
    ax.set_xlabel('機首からの距離 x [mm]' if JP_OK else 'Distance from nose x [mm]')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)

    title = (f"機体側面図 ({target_label}) | "
             f"W={total_mass:.1f}g, CG@MAC{(x_cg-wing_le_mm)/MAC_MM*100:.1f}%, "
             f"SM={sm['SM_pct']:.1f}%, Vh={tail['Vh']:.2f}, ℓh={tail['lh_mm']:.0f}mm")
    if not JP_OK:
        title = (f"Side View ({target_label}) | "
                 f"W={total_mass:.1f}g, CG@MAC{(x_cg-wing_le_mm)/MAC_MM*100:.1f}%, "
                 f"SM={sm['SM_pct']:.1f}%, Vh={tail['Vh']:.2f}, lh={tail['lh_mm']:.0f}mm")
    ax.set_title(title, fontsize=11)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  → 保存: {save_path}")


def plot_cg_sweep(save_path: str, scenario: str = "main"):
    """主翼LE位置を変化させたときのCG位置・SM・Vhのスイープ"""
    wing_le_range = np.arange(100, 180, 1)

    cg_results = []
    sm_results = []
    vh_results = []
    for wl in wing_le_range:
        components = get_components(wl, scenario)
        cg = compute_cg(components)
        sm = compute_static_margin(wl, cg)
        tail = compute_tail_volume(wl)
        cg_results.append(cg)
        sm_results.append(sm['SM_pct'])
        vh_results.append(tail['Vh'])

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    # ① CG vs wing_LE
    ax = axes[0]
    ax.plot(wing_le_range, cg_results, 'b-', lw=2,
            label='CG位置 (計算)' if JP_OK else 'CG (calculated)')
    ax.plot(wing_le_range, wing_le_range + 0.30 * MAC_MM, 'r--', alpha=0.7,
            label='30% MAC 目標' if JP_OK else 'Target 30% MAC')
    ax.plot(wing_le_range, wing_le_range + 0.35 * MAC_MM, 'g--', alpha=0.7,
            label='35% MAC 目標' if JP_OK else 'Target 35% MAC')

    # 解の交点をマーク
    wl_30 = find_wing_le_for_cg_target(0.30, scenario)
    wl_35 = find_wing_le_for_cg_target(0.35, scenario)
    cg_30 = wl_30 + 0.30 * MAC_MM
    cg_35 = wl_35 + 0.35 * MAC_MM
    ax.plot(wl_30, cg_30, 'ro', markersize=10)
    ax.plot(wl_35, cg_35, 'go', markersize=10)
    ax.annotate(f'  ({wl_30:.1f}, {cg_30:.1f})', (wl_30, cg_30), fontsize=9, color='red')
    ax.annotate(f'  ({wl_35:.1f}, {cg_35:.1f})', (wl_35, cg_35), fontsize=9, color='green')

    ax.set_ylabel('CG位置 [mm]' if JP_OK else 'CG position [mm]')
    ax.set_title('主翼LE位置 vs CG位置' if JP_OK else 'Wing LE vs CG position', fontsize=11)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # ② SM vs wing_LE
    ax = axes[1]
    ax.plot(wing_le_range, sm_results, 'purple', lw=2)
    ax.axhspan(5, 15, alpha=0.2, color='green', label='目標範囲 5-15%' if JP_OK else 'Target 5-15%')
    ax.axhline(0, color='red', linestyle='--', alpha=0.5,
               label='中立 (不安定境界)' if JP_OK else 'Neutral')
    ax.axvline(wl_30, color='red', linestyle=':', alpha=0.6)
    ax.axvline(wl_35, color='green', linestyle=':', alpha=0.6)
    ax.set_ylabel('静安定余裕 SM [%MAC]' if JP_OK else 'Static Margin [%MAC]')
    ax.set_title('主翼LE位置 vs 静安定余裕' if JP_OK else 'Wing LE vs Static Margin', fontsize=11)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # ③ Vh vs wing_LE
    ax = axes[2]
    ax.plot(wing_le_range, vh_results, 'orange', lw=2)
    ax.axhspan(0.5, 1.0, alpha=0.2, color='green',
               label='教科書範囲 0.5-1.0' if JP_OK else 'Textbook 0.5-1.0')
    ax.axhline(0.55, color='black', linestyle=':', alpha=0.7,
               label='設計値 Vh=0.55' if JP_OK else 'Design Vh=0.55')
    ax.axvline(wl_30, color='red', linestyle=':', alpha=0.6)
    ax.axvline(wl_35, color='green', linestyle=':', alpha=0.6)
    ax.set_xlabel('主翼LE位置 wing_LE [mm]' if JP_OK else 'Wing LE [mm]')
    ax.set_ylabel('容積比 Vh' if JP_OK else 'Tail volume Vh')
    ax.set_title('主翼LE位置 vs 水平尾翼容積比 Vh' if JP_OK else 'Wing LE vs Vh', fontsize=11)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  → 保存: {save_path}")


# =============================================================================
# Section 5: メイン解析
# =============================================================================

def print_separator(char='=', n=70):
    print(char * n)


def analyze_configuration(target_mac_frac: float, scenario: str = "main"):
    """1パターンの解析を実施し、結果を表示"""
    target_label = f"CG目標 MAC {int(target_mac_frac*100)}%"

    print_separator('-')
    print(f"  {target_label}, シナリオ: {WEIGHT_SCENARIOS[scenario]['label']}")
    print_separator('-')

    # 解を求める
    wing_le = find_wing_le_for_cg_target(target_mac_frac, scenario)
    components = get_components(wing_le, scenario)
    x_cg = compute_cg(components)
    total_mass = compute_total_mass(components)
    sm = compute_static_margin(wing_le, x_cg)
    tail = compute_tail_volume(wing_le)

    # 部品リスト表示
    print(f"  {'部品':<25} {'質量[g]':>8} {'x位置[mm]':>10} {'モーメント':>12}")
    for c in components:
        print(f"  {c.name:<25} {c.mass_g:>8.2f} {c.x_mm:>10.1f} {c.mass_g*c.x_mm:>12.1f}")
    print(f"  {'合計':<25} {total_mass:>8.2f}")

    # 主結果
    print()
    print(f"  ▶ 主翼LE位置:      wing_LE = {wing_le:.1f} mm")
    print(f"  ▶ 主翼TE位置:      wing_TE = {wing_le + MAC_MM:.1f} mm")
    print(f"  ▶ 重心位置:        x_cg    = {x_cg:.1f} mm  (= MAC {(x_cg-wing_le)/MAC_MM*100:.1f}%)")
    print(f"  ▶ 機体総重量:      W       = {total_mass:.1f} g")
    print(f"  ▶ 中立点:          x_np    = {sm['x_np_mm']:.1f} mm  (= MAC {X_NP_MAC_FRAC*100:.1f}%)")
    print(f"  ▶ 静安定余裕:      SM      = {sm['SM_pct']:.1f}% MAC")
    print(f"  ▶ 尾翼アーム:      ℓh      = {tail['lh_mm']:.1f} mm")
    print(f"  ▶ 水平尾翼容積比:  Vh      = {tail['Vh']:.3f}  (設計値0.55)")
    print(f"  ▶ 垂直尾翼容積比:  Vv      = {tail['Vv']:.3f}  (設計値0.04)")

    # 評価コメント
    print()
    if 5 <= sm['SM_pct'] <= 15:
        sm_eval = "✓ 目標範囲内(5-15%)"
    elif 15 < sm['SM_pct'] <= 25:
        sm_eval = "△ やや過安定(動きが鈍め、初心者向きとしては悪くない)"
    elif sm['SM_pct'] > 25:
        sm_eval = "✗ 過安定すぎ(操縦応答が鈍い)"
    else:
        sm_eval = "✗ 不安定 or 動きすぎ"
    print(f"  SM評価: {sm_eval}")

    if 0.45 <= tail['Vh'] <= 0.65:
        vh_eval = "✓ 設計値0.55の±15%以内"
    elif tail['Vh'] > 0.65:
        vh_eval = "△ 尾翼が大きめ(操舵効きすぎ気味)"
    else:
        vh_eval = "✗ 尾翼不足"
    print(f"  Vh評価: {vh_eval}")

    return {
        'wing_le': wing_le,
        'x_cg': x_cg,
        'total_mass': total_mass,
        'sm': sm['SM_pct'],
        'lh': tail['lh_mm'],
        'Vh': tail['Vh'],
        'Vv': tail['Vv'],
    }


def main():
    print_separator()
    print("  §7 重量配分・重心解析  (v3: ツインブーム+U字+ロッドフェアリング)")
    print(f"  v3更新: 2026-05-17  POWERUP実測値: {POWERUP_MASS_G}g @ x={POWERUP_CG_X_MM}mm")
    print_separator()

    # 出力ディレクトリ(このスクリプトと同じ階層の output/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # === パターン①: CG目標 MAC 30%、主シナリオ ===
    result_30 = analyze_configuration(0.30, "main")
    plot_side_view(result_30['wing_le'], "main",
                   f"{output_dir}/weight_distribution_2026-05-06_CG30.png",
                   "CG目標 MAC 30%")

    # === パターン②: CG目標 MAC 35%、主シナリオ ===
    print()
    result_35 = analyze_configuration(0.35, "main")
    plot_side_view(result_35['wing_le'], "main",
                   f"{output_dir}/weight_distribution_2026-05-06_CG35.png",
                   "CG目標 MAC 35%")

    # === CGスイーププロット ===
    print()
    print_separator('-')
    print("  CGスイーププロット生成中...")
    print_separator('-')
    plot_cg_sweep(f"{output_dir}/cg_sweep_2026-05-06.png", "main")

    # === 重量シナリオ感度解析 ===
    print()
    print_separator('-')
    print("  重量シナリオ感度解析(設計点 CG目標35%固定)")
    print_separator('-')
    print(f"  {'シナリオ':<20} {'W [g]':>8} {'wing_LE [mm]':>14} {'ℓh [mm]':>10} {'Vh':>8}")
    for scn_key in ["optimistic", "main", "pessimistic"]:
        wl = find_wing_le_for_cg_target(0.35, scn_key)
        comps = get_components(wl, scn_key)
        W = compute_total_mass(comps)
        tail = compute_tail_volume(wl)
        label = WEIGHT_SCENARIOS[scn_key]['label']
        print(f"  {label:<20} {W:>8.1f} {wl:>14.1f} {tail['lh_mm']:>10.1f} {tail['Vh']:>8.3f}")

    # === サマリ ===
    print()
    print_separator()
    print("  サマリ(主シナリオ)")
    print_separator()
    print(f"  ▶ CG@30% MAC: wing_LE={result_30['wing_le']:.1f}mm, "
          f"W={result_30['total_mass']:.1f}g, SM={result_30['sm']:.1f}%, Vh={result_30['Vh']:.3f}")
    print(f"  ▶ CG@35% MAC: wing_LE={result_35['wing_le']:.1f}mm, "
          f"W={result_35['total_mass']:.1f}g, SM={result_35['sm']:.1f}%, Vh={result_35['Vh']:.3f}")
    print()
    print(f"  生成ファイル({output_dir}):")
    print(f"    - weight_distribution_2026-05-06_CG30.png")
    print(f"    - weight_distribution_2026-05-06_CG35.png")
    print(f"    - cg_sweep_2026-05-06.png")


if __name__ == "__main__":
    main()