"""
mainwing_structural_analysis.py
================================

§2 主翼の概略設計 — 空力・構造の数値計算スクリプト

対象 : 5/18 第2回審査会の概略計算書 §2
作成 : 2026-05-05
v3   : 2026-05-17 — 機体諸元(W=36.5 g、巡航 V=7.5 m/s、|L_t|=39.8 mN)を
       §3 尾翼U字版v3・§7 重量重心v3 と統一して再整理。

役割 : §2 主翼計算書の各数値の「根拠プログラム」。
       XFLR5 v6.62 の VLM 3次元解析(α=3°, V=7 m/s, AG36)の生データを読み、
       巡航運用点(V=7.5 m/s)へ外挿して翼根曲げモーメント・応力・安全率を出す。

入力 : vlm_data/spanwise_alpha3deg_AR6_V7ms_2026-05-05.txt(XFLR5 OpPoint Export)
出力 : 標準出力に VLM 生データ・巡航外挿・翼根応力・安全率

注記 : VLM は α=3°, V=7 m/s の1点で実施済み。巡航運用点(α≈3.7°, V=7.5 m/s)の
       値はこの基準点からの線形外挿であり、§2.7 の残課題で運用点直接解析を予定。
"""

import os
import re

# =============================================================================
# 入力定数
# =============================================================================

HERE = os.path.dirname(os.path.abspath(__file__))
VLM_FILE = os.path.join(HERE, "vlm_data", "spanwise_alpha3deg_AR6_V7ms_2026-05-05.txt")

RHO = 1.225          # 空気密度 [kg/m³]
G = 9.81             # 重力加速度 [m/s²]

# --- 主翼諸元(§2)---
S_WING = 216e-4      # 翼面積 [m²]
B_WING = 0.360       # 翼幅 [m]
C_WING = 0.060       # 翼弦 [m]

# --- 機体諸元(§3・§7 U字版v3 と統一)---
W_TOTAL = 0.358      # 機体総重量 [N](36.5 g)
L_TAIL = 0.0398      # 水平尾翼トリム荷重 |L_t| [N](§3 v3:下向き 39.8 mN)
V_CRUISE = 7.5       # 巡航速度 [m/s]

LOAD_FACTOR = 2.0    # 荷重倍率 n(突風・操舵、§2.1)

# --- スパー(主桁)---
SPAR_DIAMETER = 0.002      # φ2 mm CFRP ロッド [m]
SIGMA_ULT_CFRP = 1500e6    # CFRP 引張強度 [Pa](代表値、§2.7 で実購入品に更新)


def dynamic_pressure(v, rho=RHO):
    return 0.5 * rho * v ** 2


def section_modulus_circular(d):
    """円形断面の断面係数 W = π·d³/32"""
    import math
    return math.pi * d ** 3 / 32.0


def parse_vlm_header(path):
    """XFLR5 OpPoint Export のヘッダから CL/Cd/Cm/Bending 等を抽出する。"""
    vals = {}
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            for key in ('QInf', 'Alpha', 'CL', 'Cd', 'ICd', 'PCd', 'Cm', 'Bending'):
                m = re.search(rf'(?<![A-Za-z]){key}\s*=\s*([-\d.eE+]+)', line)
                if m and key not in vals:
                    vals[key] = float(m.group(1))
            if 'Main Wing' in line:
                break
    return vals


def main():
    print("=" * 72)
    print(" §2 主翼の概略設計  空力・構造の数値計算")
    print("=" * 72)

    if not os.path.exists(VLM_FILE):
        print(f"[エラー] VLM データ未発見: {VLM_FILE}")
        return

    v = parse_vlm_header(VLM_FILE)
    q_vlm = dynamic_pressure(v['QInf'])
    LD = v['CL'] / v['Cd']

    print(f"\n[VLM 3次元解析 生データ]  XFLR5 v6.62 / AG36 / AR=6")
    print(f"  ファイル: {os.path.basename(VLM_FILE)}")
    print(f"  解析条件 : α = {v['Alpha']:.1f}°,  V = {v['QInf']:.1f} m/s,  q = {q_vlm:.2f} Pa")
    print(f"  CL       = {v['CL']:.4f}")
    print(f"  CD       = {v['Cd']:.4f}  (誘導 ICd={v['ICd']:.4f} + 粘性 PCd={v['PCd']:.4f})")
    print(f"  Cm       = {v['Cm']:.5f}   (水平尾翼設計の入力値)")
    print(f"  L/D      = {LD:.2f}")
    print(f"  翼根曲げモーメント Bending = {v['Bending'] * 1000:.3f} mN·m  (XFLR5 自動算出)")

    # VLM 基準点の半翼揚力
    L_vlm = v['CL'] * q_vlm * S_WING

    # ---------- 巡航運用点への外挿 ----------
    print(f"\n[巡航運用点への外挿]  V = {V_CRUISE} m/s")
    q_cruise = dynamic_pressure(V_CRUISE)
    L_wing_cruise = W_TOTAL + L_TAIL
    CL_cruise = L_wing_cruise / (q_cruise * S_WING)
    print(f"  動圧 q              = {q_cruise:.2f} Pa")
    print(f"  必要主翼揚力 L_wing = W + |L_t| = {W_TOTAL:.3f} + {L_TAIL:.4f} = {L_wing_cruise:.4f} N")
    print(f"  必要揚力係数 CL     = L_wing/(q·S) = {CL_cruise:.3f}")

    # ---------- 翼根曲げモーメント ----------
    # M_root は揚力の大きさに比例する(スパン荷重形状は α でほぼ不変)。
    # VLM 基準点の Bending を、揚力比でスケールして巡航運用点の値を求める。
    lift_ratio = L_wing_cruise / L_vlm
    M_root_n1 = v['Bending'] * lift_ratio
    M_root_n2 = M_root_n1 * LOAD_FACTOR

    print(f"\n[翼根曲げモーメント]  M_root = Bending_VLM × 揚力比")
    print(f"  VLM 基準点の主翼揚力 L_vlm = CL·q·S = {L_vlm:.4f} N")
    print(f"  揚力比 = L_wing / L_vlm    = {lift_ratio:.4f}")
    print(f"  巡航時 (n=1)  M_root = {M_root_n1 * 1000:.2f} mN·m")
    print(f"  設計時 (n={LOAD_FACTOR:.0f})  M_root = {M_root_n2 * 1000:.2f} mN·m")

    # ---------- 翼根応力・安全率 ----------
    Ws = section_modulus_circular(SPAR_DIAMETER)
    sigma = M_root_n2 / Ws
    SF = SIGMA_ULT_CFRP / sigma

    print(f"\n[翼根応力・安全率]  φ{SPAR_DIAMETER*1000:.0f} mm CFRP スパー")
    print(f"  断面係数 W_s = π·d³/32 = {Ws * 1e9:.3f} mm³")
    print(f"  翼根応力 σ   = M_root(n=2) / W_s = {sigma / 1e6:.1f} MPa")
    print(f"  安全率   SF  = σ_ult / σ = {SIGMA_ULT_CFRP/1e6:.0f} / {sigma/1e6:.1f} = {SF:.0f}")

    # ---------- サマリー ----------
    print("\n" + "=" * 72)
    print(" 確定値サマリー(計算書 §2 提供物)")
    print("=" * 72)
    print(f"  VLM 基準(α=3°,V=7): CL={v['CL']:.4f}  CD={v['Cd']:.4f}  "
          f"L/D={LD:.2f}  Cm={v['Cm']:.4f}")
    print(f"  巡航運用点(V=7.5) : CL={CL_cruise:.3f}")
    print(f"  翼根曲げ M_root    : 巡航 {M_root_n1*1000:.2f} / 設計 {M_root_n2*1000:.2f} mN·m")
    print(f"  翼根応力・安全率   : σ={sigma/1e6:.1f} MPa,  SF={SF:.0f}")
    print("\n[完了]")


if __name__ == '__main__':
    main()
