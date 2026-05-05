"""
主翼構造解析スクリプト
=======================

XFLR5 VLM 3D解析の翼幅方向データを読み込んで、以下を一気に実施:
  ① Cl分布グラフ作成(計算書貼付用)
  ② 翼根曲げモーメントの手計算(XFLR5検算)
  ③ カーボンスパー翼根曲げ応力計算
  ④ 安全率評価

入力:
  vlm_data/spanwise_alpha3deg_AR6_V7ms_2026-05-05.txt

出力:
  plots/cl_distribution_alpha3deg_cruise.png
  plots/bending_moment_alpha3deg.png
  コンソール出力: 主要数値一式(計算書に貼り付ける用)

実行:
  python mainwing_structural_analysis.py

作成: 2026-05-05
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 日本語フォント設定(Windows標準フォントを優先順に指定)
rcParams['font.family'] = ['Meiryo', 'Yu Gothic', 'MS Gothic', 'sans-serif']
rcParams['axes.unicode_minus'] = False  # マイナス記号の文字化け防止

# numpy 2.x では np.trapz が削除されたため、np.trapezoid に統一
# (古いバージョンでも動くようフォールバック)
try:
    _trapz = np.trapezoid
except AttributeError:
    _trapz = np.trapz

# ============================================================
# 設定値
# ============================================================

# パス
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(HERE, 'vlm_data',
                         'spanwise_alpha3deg_AR6_V7ms_2026-05-05.txt')
PLOTS_DIR = os.path.join(HERE, 'plots')
os.makedirs(PLOTS_DIR, exist_ok=True)

# 飛行条件 (要件定義書 §4.2-4.3)
V = 7.0  # 巡航速度 [m/s]
RHO = 1.225  # 空気密度 [kg/m³]
S = 0.0216  # 翼面積 [m²]
B = 0.360  # 翼幅 [m]
CHORD = 0.060  # 翼弦 [m] (矩形なので一定)
MASS = 0.0325  # 機体総重量 [kg]
G = 9.81  # 重力加速度 [m/s²]

# カーボンスパー諸元 (φ2mm 円柱断面 pultruded CFRP)
SPAR_D = 0.002  # 直径 [m]
SPAR_I = np.pi * SPAR_D ** 4 / 64  # 断面二次モーメント [m⁴]
SPAR_C = SPAR_D / 2  # 中立軸からの最大距離 [m]
SPAR_W = SPAR_I / SPAR_C  # 断面係数 [m³]

# CFRP 材料特性(典型的な pultruded carbon rod)
# 注:正式な材料データシート入手後に要更新
SIGMA_ULT = 1500e6  # 引張強度 [Pa] = 1500 MPa
SF_DESIGN = 2.0  # 設計安全率
SIGMA_ALLOW = SIGMA_ULT / SF_DESIGN  # 許容応力 [Pa] = 750 MPa

# 飛行荷重倍数(突風・旋回時の最大加速度倍率)
LOAD_FACTOR = 2.0


# ============================================================
# データ読み込み
# ============================================================

def load_xflr5_spanwise(filepath):
    """
    XFLR5 OpPoint Export形式の.txtから、ヘッダー値と翼幅方向データ表を抽出。

    Args:
        filepath: .txtファイルパス

    Returns:
        df: DataFrame(38行 × 12列のスパン方向データ)
        header: dict(V, alpha, CL, Cm, M_root_xflr5 等)
    """
    with open(filepath, encoding='utf-8') as f:
        lines = f.readlines()

    # ヘッダー部から主要数値を抽出
    header = {}
    for line in lines[:30]:
        s = line.strip()
        if s.startswith('QInf'):
            header['V'] = float(s.split('=')[1].split()[0])
        elif s.startswith('Alpha'):
            header['alpha'] = float(s.split('=')[1].split()[0])
        elif s.startswith('CL '):
            header['CL'] = float(s.split('=')[1].strip())
        elif s.startswith('Cd '):
            # "Cd = 0.045414  ICd = ..." の形式
            header['CD'] = float(s.split('=')[1].split('ICd')[0].strip())
        elif s.startswith('Cm '):
            header['Cm'] = float(s.split('=')[1].strip())
        elif s.startswith('Bending'):
            header['M_root_xflr5'] = float(s.split('=')[1].strip())

    # 翼幅方向データ表の位置特定(ヘッダー行を探す)
    header_idx = None
    for i, line in enumerate(lines):
        if 'y-span' in line and 'Chord' in line and 'BM' in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("スパン方向データのヘッダーが見つかりません")

    # データ行を抽出(空行 or 'Main Wing Cp' で終了)
    data_rows = []
    for line in lines[header_idx + 1:]:
        s = line.strip()
        if not s or 'Cp' in s:
            break
        parts = s.split()
        if len(parts) >= 12:
            data_rows.append([float(x) for x in parts])

    columns = ['y_span', 'chord', 'Ai', 'Cl', 'PCd', 'ICd',
               'CmGeom', 'CmAirf', 'XTrtop', 'XTrBot', 'XCP', 'BM']
    df = pd.DataFrame(data_rows, columns=columns)
    return df, header


# ============================================================
# プロット作成
# ============================================================

def plot_cl_distribution(df, header, save_path):
    """Cl分布プロット(計算書用、見やすいスタイル)"""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(df['y_span'] * 1000, df['Cl'], 'b-o', markersize=4,
            linewidth=1.5, label='AG36 主翼')
    ax.axhline(header['CL'], color='red', linestyle='--', alpha=0.6,
               label=f'平均 $C_L$ = {header["CL"]:.3f}')
    ax.axvline(0, color='gray', linestyle=':', alpha=0.4)

    ax.set_xlabel('翼幅方向位置 $y$ [mm]', fontsize=11)
    ax.set_ylabel('局所揚力係数 $C_l$', fontsize=11)
    ax.set_title(f'翼幅方向 $C_l$ 分布 ($\\alpha$=3°、巡航時)\n'
                 f'$V$={V} m/s、AR=6、AG36、Re$\\approx$28,000', fontsize=11)
    ax.set_ylim(0, max(df['Cl']) * 1.15)
    ax.set_xlim(-200, 200)
    ax.legend(loc='lower center', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_bending_moment(y, M_manual, M_xflr5, save_path):
    """曲げモーメント分布プロット(検算結果込み)"""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(y * 1000, M_manual * 1000, 'g-o', markersize=4,
            linewidth=1.5, label='手計算による積分(本スクリプト)')
    ax.plot(y * 1000, M_xflr5 * 1000, 'r--', linewidth=1.5,
            alpha=0.7, label='XFLR5 BM 出力')

    ax.set_xlabel('翼幅方向位置 $y$ [mm]', fontsize=11)
    ax.set_ylabel('曲げモーメント $M(y)$ [mN·m]', fontsize=11)
    ax.set_title(f'翼幅方向 曲げモーメント分布 ($\\alpha$=3°)\n'
                 f'半翼、翼根 $y$=0', fontsize=11)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================
# メイン処理
# ============================================================

def main():
    print("=" * 60)
    print(" 主翼構造解析(α=3° 巡航条件)")
    print("=" * 60)

    # データ読み込み
    df, header = load_xflr5_spanwise(DATA_FILE)

    print("\n--- XFLR5解析結果(ヘッダーから抽出) ---")
    print(f"  速度 V        = {header['V']:.2f} m/s")
    print(f"  迎角 α        = {header['alpha']:.2f}°")
    print(f"  CL           = {header['CL']:.4f}")
    print(f"  CD           = {header['CD']:.4f}")
    print(f"  Cm           = {header['Cm']:.4f}")
    print(f"  L/D          = {header['CL'] / header['CD']:.2f}")
    print(f"  M_root(XFLR5)= {header['M_root_xflr5'] * 1000:.3f} mN·m")
    print(f"  スパン方向データ点数: {len(df)}")

    # 動圧と必要CL確認
    q = 0.5 * RHO * V ** 2
    W = MASS * G
    CL_required = 2 * W / (RHO * V ** 2 * S)
    print(f"\n--- 飛行条件チェック ---")
    print(f"  動圧 q = ½ρV²    = {q:.2f} Pa")
    print(f"  機体重量 W       = {W:.4f} N ({MASS * 1000:.1f} g)")
    print(f"  必要CL(水平飛行) = {CL_required:.4f}")
    print(f"  実際のCL(XFLR5)  = {header['CL']:.4f}")
    print(
        f"  → 釣り合い誤差: {abs(CL_required - header['CL']) / CL_required * 100:.1f}%")

    # ============================================================
    # ① Cl分布プロット
    # ============================================================
    print("\n--- ① Cl分布プロット ---")
    plot1_path = os.path.join(PLOTS_DIR, 'cl_distribution_alpha3deg_cruise.png')
    plot_cl_distribution(df, header, plot1_path)
    print(f"  保存: {plot1_path}")

    Cl_max = df['Cl'].max()
    Cl_root = df.iloc[len(df) // 2]['Cl']
    Cl_tip_outer = df['Cl'].iloc[-1]
    print(f"  Cl 最大値(中央付近): {Cl_max:.4f}")
    print(f"  Cl 翼端値          : {Cl_tip_outer:.4f}")
    print(f"  → 矩形翼AR=6の典型的な凸型分布")

    # ============================================================
    # ② 曲げモーメント手計算(検算)
    # ============================================================
    print("\n--- ② 曲げモーメント手計算 ---")

    # 半翼に着目: y >= 0 のデータのみ
    df_half = df[df['y_span'] >= 0].sort_values('y_span').reset_index(drop=True)
    y = df_half['y_span'].values
    Cl = df_half['Cl'].values
    chord = df_half['chord'].values

    # 局所揚力(単位スパンあたり) [N/m]
    dL_dy = q * chord * Cl

    # 半翼揚力(検算)
    L_half = _trapz(dL_dy, y)

    # 各位置 y0 での曲げモーメント:
    #   M(y0) = ∫_{y0}^{b/2} q × c(y) × Cl(y) × (y - y0) dy
    # = 翼端方向の揚力が y0 に対して作るモーメント
    M_distribution = np.zeros_like(y)
    for i in range(len(y)):
        y_outer = y[i:]
        dL_outer = dL_dy[i:]
        arms = y_outer - y[i]
        M_distribution[i] = _trapz(dL_outer * arms, y_outer)

    M_root_manual = M_distribution[0]

    print(
        f"  半翼揚力 L_half     = {L_half:.4f} N (期待値 W/2 = {W / 2:.4f} N)")
    print(f"  → 釣り合い誤差     : {abs(L_half - W / 2) / (W / 2) * 100:.1f}%")
    print(f"  翼根曲げモーメント:")
    print(f"    手計算 M_root     = {M_root_manual * 1000:.3f} mN·m")
    print(f"    XFLR5  M_root     = {header['M_root_xflr5'] * 1000:.3f} mN·m")
    diff_pct = abs(M_root_manual - header['M_root_xflr5']) / header[
        'M_root_xflr5'] * 100
    print(f"    差               : {diff_pct:.2f}%")
    if diff_pct < 5.0:
        print(f"    → ✓ 検算OK(差 < 5%)")
    else:
        print(f"    → ⚠ 差が大きい、要確認")

    # 曲げモーメント分布プロット
    plot2_path = os.path.join(PLOTS_DIR, 'bending_moment_alpha3deg.png')
    plot_bending_moment(y, M_distribution, df_half['BM'].values, plot2_path)
    print(f"  保存: {plot2_path}")

    # ============================================================
    # ③ 翼根曲げ応力計算(φ2mm カーボンスパー)
    # ============================================================
    print("\n--- ③ 翼根曲げ応力計算(φ2mm CFRPスパー) ---")
    print(f"  スパー諸元:")
    print(f"    直径 d        = {SPAR_D * 1000:.1f} mm")
    print(f"    断面二次モーメント I = πd⁴/64 = {SPAR_I * 1e12:.4f} mm⁴")
    print(f"    断面係数 W = I/c     = {SPAR_W * 1e9:.4f} mm³")

    # 巡航時(荷重倍率1)の曲げ応力
    M_design_cruise = M_root_manual  # 巡航時
    M_design_max = M_root_manual * LOAD_FACTOR  # 設計時(突風・旋回想定)

    sigma_cruise = M_design_cruise / SPAR_W
    sigma_design = M_design_max / SPAR_W

    print(f"\n  巡航時(荷重倍率 n=1):")
    print(f"    σ_cruise  = M/W = {sigma_cruise / 1e6:.2f} MPa")
    print(f"  設計時(荷重倍率 n={LOAD_FACTOR}):")
    print(f"    σ_design  = M×n/W = {sigma_design / 1e6:.2f} MPa")

    # ============================================================
    # ④ 安全率評価
    # ============================================================
    print("\n--- ④ 安全率評価 ---")
    print(f"  CFRP材料特性(典型値、要更新):")
    print(f"    引張強度 σ_ult     = {SIGMA_ULT / 1e6:.0f} MPa")
    print(
        f"    設計許容応力 σ_allow = σ_ult/{SF_DESIGN} = {SIGMA_ALLOW / 1e6:.0f} MPa")

    SF_cruise = SIGMA_ULT / sigma_cruise
    SF_design = SIGMA_ULT / sigma_design

    print(f"\n  安全率(対 σ_ult):")
    print(f"    巡航時 SF       = {SF_cruise:.1f}")
    print(f"    設計時 SF(n={LOAD_FACTOR}) = {SF_design:.1f}")

    print(f"\n  判定:")
    if sigma_design < SIGMA_ALLOW:
        margin_pct = (SIGMA_ALLOW - sigma_design) / SIGMA_ALLOW * 100
        print(
            f"    σ_design ({sigma_design / 1e6:.2f} MPa) < σ_allow ({SIGMA_ALLOW / 1e6:.0f} MPa)")
        print(f"    余裕度: {margin_pct:.1f}%")
        print(f"    → ✓ 設計成立")
        if SF_design > 20:
            print(f"    補足: SF={SF_design:.0f} は過剰の可能性。")
            print(f"          スパー直径削減(φ1.5mm 等)で軽量化可能性あり")
    else:
        print(
            f"    σ_design ({sigma_design / 1e6:.2f} MPa) >= σ_allow ({SIGMA_ALLOW / 1e6:.0f} MPa)")
        print(f"    → ✗ 設計不成立 — スパー強化 or 材質変更要")

    print("\n" + "=" * 60)
    print(" 解析完了")
    print("=" * 60)


if __name__ == '__main__':
    main()