"""
Tail_design_analysis.py
========================

§3 尾翼の概略設計 — 数値計算スクリプト(v3:ツインブーム + U字型尾翼)

対象 : 5/18 第2回審査会の概略計算書 §3
作成 : 2026-05-07
v3   : 2026-05-17 — ツインブーム + U字型尾翼 + ロッドフェアリング構成へ更新。
       §2 主翼・§5 プロペラ・§7 重量重心 U字版v3 と数値を統一。

役割 : §3 尾翼U字版v3 計算書の各数値の「根拠プログラム」。
入力 : 機体諸元(§2 主翼・§7 重量重心v3)、xfoil_data/NACA0008 ポーラ
出力 : 標準出力に Step ①〜③ の計算結果

カバー範囲:
  Step ①   力・モーメントの設定(トリム尾翼力 L_t、最悪荷重)
  Step ②-1 容積比からの尾翼面積 Sh, Sv
  Step ②-2 寸法(水平尾翼=ブーム間隔固定、垂直尾翼=U字2枚)
  Step ②-3 翼型 NACA0008 ポーラ照合
  Step ③   強度設計(水平尾翼・垂直尾翼・ツインブーム)
  + 静安定余裕 SM
"""

import os
import numpy as np


# =============================================================================
# 入力定数(全て §2 主翼 / §7 重量重心 U字版v3 由来)
# =============================================================================

# --- 環境定数 ---
RHO = 1.225          # 空気密度 [kg/m³](海面標準)
NU = 1.5e-5          # 動粘性係数 [m²/s]
G = 9.81             # 重力加速度 [m/s²]

# --- 主翼仕様(§2)---
S_WING = 216e-4      # 主翼面積 [m²](216 cm²)
B_WING = 0.360       # 主翼幅 [m]
C_WING = 0.060       # MAC = 翼弦 [m]
AR_WING = 6.0        # 主翼アスペクト比

# --- 巡航条件(§2、VLM α=3°/V=7 を巡航点 V=7.5 へ外挿)---
V_CRUISE = 7.5       # 巡航速度 [m/s]
ALPHA_CRUISE = 3.7   # 巡航迎角 [deg](VLM α=3°基準からの外挿)
CL_CRUISE = 0.535    # 巡航 CL(L_wing=0.398 N @ V=7.5)
CM_WING = -0.162     # 主翼ピッチングモーメント係数(VLM -0.16238 を3桁丸め、主翼AC基準)

# --- 機体配置(§7 重量重心 U字版v3 で確定)---
W_MAIN = 36.47e-3    # 機体総重量 [kg] 主シナリオ

# 設計点(CG @ MAC 35%)
WING_LE_DESIGN = 0.1280   # 主翼前縁位置 [m]
X_AC_DESIGN    = WING_LE_DESIGN + 0.25 * C_WING   # 主翼AC = 0.1430 m
X_CG_DESIGN    = 0.14894  # 機体CG位置 [m](§7 v3)

# 限界点(CG @ MAC 30%、重量バラつき余裕の許容端)
WING_LE_LIMIT  = 0.1315
X_AC_LIMIT     = WING_LE_LIMIT + 0.25 * C_WING    # = 0.1465 m
X_CG_LIMIT     = 0.14942  # 機体CG位置 [m](§7 v3、主翼後方スライド後)

# 尾翼AC位置 — v3 機体配置のアンカー(機体最後端 305 mm、c_h≈37→36 でも不変)
X_TAIL_AC      = 0.27725  # 水平尾翼AC位置 [m]
ELL_H_DESIGN   = X_TAIL_AC - X_AC_DESIGN   # 主翼AC↔尾翼AC 設計点 [m] = 0.13425
ELL_H_LIMIT    = X_TAIL_AC - X_AC_LIMIT    # 限界点 [m] = 0.13075

# --- 容積比目標(Drela MIT 16.01 "Lab 8 Notes")---
VH_TARGET = 0.55     # 水平尾翼容積比(推奨 0.30〜0.60 の上端寄り、屋内+初心者で安定側)
VV_TARGET = 0.04     # 垂直尾翼容積比(推奨 0.02〜0.05 の中央〜上端寄り)

# --- 寸法の確定値(v3)---
B_BOOM = 0.148       # ツインブーム間隔 = 水平尾翼翼幅 [m]
C_TAIL_H = 0.036     # 水平尾翼翼弦 [m](Sh/B_BOOM を切れの良い値へ丸め)
B_VTAIL = 0.040      # 垂直尾翼 翼幅(高さ)[m] /枚
C_VTAIL = 0.030      # 垂直尾翼 翼弦 [m] /枚
N_VTAIL = 2          # U字型 垂直尾翼 枚数
ELL_BOOM = 0.200     # ツインブーム全長 [m](保守側、片持ち長)

# --- 安全係数 ---
LOAD_FACTOR = 2.0    # 荷重倍率 n(突風・操舵、§2 主翼と統一)
L_WORST_DESIGN = 0.100   # 強度設計荷重 [N](実限界 96.8 mN を保守側に据え置き)

# --- 材料定数 ---
SIGMA_ULT_CFRP  = 1500e6   # CFRP 引張強度 [Pa](代表値)
SIGMA_YIELD_EPS = 0.4e6    # EPS 降伏応力 [Pa](低密度フォーム代表値)


# =============================================================================
# 関数定義
# =============================================================================

def dynamic_pressure(v, rho=RHO):
    """動圧 q = ½ρV² [Pa]"""
    return 0.5 * rho * v ** 2


def reynolds_number(v, c, nu=NU):
    """レイノルズ数 Re = Vc/ν"""
    return v * c / nu


def trim_tail_force(cm_wing, q, S, c, W_force, x_cg, x_ac, x_tail_ac):
    """
    機体CGまわりのモーメント釣り合いから巡航時の尾翼力 L_t を計算。

      ΣM_cg = M_w,AC + W·(x_cg − x_AC) − L_t·(x_tail_ac − x_cg) = 0
      → L_t = [M_w,AC + W·(x_cg − x_AC)] / (x_tail_ac − x_cg)

    符号:L_t < 0 → 下向き(ダウンフォース)
    """
    M_ac = cm_wing * q * S * c
    M_lever = W_force * (x_cg - x_ac)
    M_net = M_ac + M_lever
    arm = x_tail_ac - x_cg          # = ℓ_cg(重心〜尾翼ACのモーメントアーム)
    L_t = M_net / arm
    return L_t, M_ac, M_lever, arm


def tail_area_from_volume(V_target, S_wing, ref_length, ell):
    """容積比から尾翼面積を逆算。 S = V·S_wing·ref / ℓ"""
    return V_target * S_wing * ref_length / ell


def neutral_point_frac(vh, ar=AR_WING, ar_h=None):
    """中立点位置 x_np/MAC(Drela MIT 16.01 "Lab 8 Notes" 近似式)"""
    return 0.25 + (1.0 - 4.0 / (ar + 2.0)) * vh * (1.0 + 2.0 / ar) / (1.0 + 2.0 / ar_h)


def section_modulus_circular(diameter):
    """円形断面の断面係数 W = π·d³/32"""
    return np.pi * diameter ** 3 / 32.0


def section_modulus_rectangle(width, height):
    """矩形断面の断面係数 W = b·h²/6"""
    return width * height ** 2 / 6.0


def root_moment_center_supported(load, span):
    """中央支持・両側張り出し梁(水平尾翼)の翼根曲げモーメント M = L·b/8"""
    return load * span / 8.0


def load_xfoil_polar(filepath):
    """XFLR5 / XFoil の export .txt をパースして dict を返す。"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    meta = {}
    for line in lines:
        if 'Re =' in line and 'Ncrit' in line:
            tokens = line.replace('=', ' ').split()
            for i, tok in enumerate(tokens):
                if tok == 'Re':
                    val = float(tokens[i + 1])
                    if i + 3 < len(tokens) and tokens[i + 2] == 'e':
                        val *= 10 ** float(tokens[i + 3])
                    meta['Re'] = val
            break
    data_rows = []
    for line in lines:
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            float(parts[0])
            data_rows.append([float(p) for p in parts[:5]])
        except ValueError:
            continue
    data = np.array(data_rows)
    return {'meta': meta, 'alpha': data[:, 0], 'CL': data[:, 1],
            'CD': data[:, 2], 'Cm': data[:, 4], 'filepath': filepath}


# =============================================================================
# メイン計算
# =============================================================================

def main():
    print("=" * 72)
    print(" §3 尾翼の概略設計  数値計算  (v3: ツインブーム + U字型尾翼)")
    print("=" * 72)

    q = dynamic_pressure(V_CRUISE)
    W_force = W_MAIN * G

    print(f"\n[環境・基本条件]")
    print(f"  巡航速度 V    = {V_CRUISE} m/s")
    print(f"  動圧     q    = {q:.2f} Pa")
    print(f"  機体重力 W    = {W_force * 1000:.1f} mN ({W_MAIN * 1000:.2f} g)")
    print(f"  主翼 Cm       = {CM_WING}")

    # ---------- Step ① 力・モーメントの設定 ----------
    print("\n" + "─" * 72)
    print(" Step ①  力・モーメントの設定")
    print("─" * 72)

    L_t_d, M_ac_d, M_lever_d, arm_d = trim_tail_force(
        CM_WING, q, S_WING, C_WING, W_force, X_CG_DESIGN, X_AC_DESIGN, X_TAIL_AC)
    print(f"\n  [設計点 CG@MAC35%]")
    print(f"    M_w,AC = Cm·q·S·c       = {M_ac_d * 1000:+7.2f} mN·m")
    print(f"    W·(x_cg−x_AC)           = {M_lever_d * 1000:+7.2f} mN·m")
    print(f"    残差(打ち消すべき M)    = {(M_ac_d + M_lever_d) * 1000:+7.2f} mN·m")
    print(f"    ℓ_cg(重心〜尾翼AC)     = {arm_d * 1000:7.2f} mm")
    print(f"    巡航尾翼力 L_t          = {L_t_d * 1000:+7.2f} mN  ({L_t_d / G * 1000:+5.2f} g重)")

    L_t_l, M_ac_l, M_lever_l, arm_l = trim_tail_force(
        CM_WING, q, S_WING, C_WING, W_force, X_CG_LIMIT, X_AC_LIMIT, X_TAIL_AC)
    print(f"\n  [限界点 CG@MAC30%]")
    print(f"    M_w,AC                  = {M_ac_l * 1000:+7.2f} mN·m")
    print(f"    W·(x_cg−x_AC)           = {M_lever_l * 1000:+7.2f} mN·m")
    print(f"    残差                    = {(M_ac_l + M_lever_l) * 1000:+7.2f} mN·m")
    print(f"    ℓ_cg                    = {arm_l * 1000:7.2f} mm")
    print(f"    巡航尾翼力 L_t          = {L_t_l * 1000:+7.2f} mN  ({L_t_l / G * 1000:+5.2f} g重)")

    L_t_worst = abs(L_t_l) * LOAD_FACTOR
    print(f"\n  [強度設計用 最悪荷重]")
    print(f"    |L_t|_worst(実限界) = 限界点 × n={LOAD_FACTOR:.0f} = {L_t_worst * 1000:.2f} mN")
    print(f"    設計荷重 L_worst(保守値)            = {L_WORST_DESIGN * 1000:.0f} mN(据え置き)")

    # ---------- Step ②-1 容積比から尾翼面積 ----------
    print("\n" + "─" * 72)
    print(" Step ②-1  容積比からの尾翼面積")
    print("─" * 72)

    Sh = tail_area_from_volume(VH_TARGET, S_WING, C_WING, ELL_H_DESIGN)
    Sv = tail_area_from_volume(VV_TARGET, S_WING, B_WING, ELL_H_DESIGN)
    print(f"\n  [設計点 ℓ_h = ℓ_v = {ELL_H_DESIGN * 1000:.2f} mm]")
    print(f"    Sh = Vh·S·c/ℓ_h = {VH_TARGET}×216×60/{ELL_H_DESIGN*1000:.2f} = {Sh * 1e4:.2f} cm²")
    print(f"    Sv = Vv·S·b/ℓ_v = {VV_TARGET}×216×360/{ELL_H_DESIGN*1000:.2f} = {Sv * 1e4:.2f} cm²")

    # ---------- Step ②-2 寸法 ----------
    print("\n" + "─" * 72)
    print(" Step ②-2  寸法(水平尾翼=ブーム間隔固定、垂直尾翼=U字2枚)")
    print("─" * 72)

    c_h_calc = Sh / B_BOOM
    Sh_adopt = B_BOOM * C_TAIL_H
    AR_h = B_BOOM / C_TAIL_H
    print(f"\n  [水平尾翼 1枚、左右ブーム間に橋渡し]")
    print(f"    翼幅 b_h = ブーム間隔 = {B_BOOM * 1000:.0f} mm(固定)")
    print(f"    必要翼弦 c_h = Sh/b_h = {Sh * 1e4:.2f}/{B_BOOM*100:.1f} = {c_h_calc * 1000:.1f} mm")
    print(f"    採用 c_h = {C_TAIL_H * 1000:.0f} mm  → 採用 Sh = {Sh_adopt * 1e4:.2f} cm²,  AR_h = {AR_h:.3f}")
    print(f"    動作 Re = {reynolds_number(V_CRUISE, C_TAIL_H):.0f}")

    Sv_adopt = N_VTAIL * B_VTAIL * C_VTAIL
    AR_v = B_VTAIL / C_VTAIL
    print(f"\n  [垂直尾翼 U字型 {N_VTAIL}枚、左右ブーム後端から上向き]")
    print(f"    1枚 {B_VTAIL*1000:.0f}×{C_VTAIL*1000:.0f} mm  → 合計 Sv = {Sv_adopt * 1e4:.1f} cm²,  AR_v = {AR_v:.3f}/枚")
    print(f"    動作 Re = {reynolds_number(V_CRUISE, C_VTAIL):.0f}")

    Vh_real = Sh_adopt * ELL_H_DESIGN / (S_WING * C_WING)
    Vv_real = Sv_adopt * ELL_H_DESIGN / (S_WING * B_WING)
    Vh_real_lim = Sh_adopt * ELL_H_LIMIT / (S_WING * C_WING)
    print(f"\n  [実現容積比(採用寸法で逆算)]")
    print(f"    Vh = Sh·ℓ_h/(S·c) = {Vh_real:.3f}(設計点) / {Vh_real_lim:.3f}(限界点)  目標 {VH_TARGET}")
    print(f"    Vv = Sv·ℓ_v/(S·b) = {Vv_real:.4f}(設計点)  目標 {VV_TARGET}")

    # ---------- Step ②-3 翼型ポーラ照合 ----------
    print("\n" + "─" * 72)
    print(" Step ②-3  翼型 NACA0008 ポーラ照合")
    print("─" * 72)

    CL_tail_d = L_t_d / (q * Sh_adopt)
    CL_tail_l = L_t_l / (q * Sh_adopt)
    print(f"\n  必要尾翼揚力係数 CL_tail = L_t/(q·Sh)")
    print(f"    設計点: {CL_tail_d:+.3f}    限界点: {CL_tail_l:+.3f}")

    polar_dir = 'xfoil_data'
    primary_name = 'NACA0008_Re16k_M0.00_N9.0.txt'
    primary_path = os.path.join(polar_dir, primary_name)
    if os.path.exists(primary_path):
        p = load_xfoil_polar(primary_path)
        mask = (p['alpha'] >= -4.5) & (p['alpha'] <= 4.5)
        a_lin, cl_lin = p['alpha'][mask], p['CL'][mask]
        slope = np.polyfit(a_lin, cl_lin, 1)[0]
        alpha_d = float(np.interp(CL_tail_d, cl_lin[np.argsort(cl_lin)],
                                  a_lin[np.argsort(cl_lin)]))
        print(f"  主シナリオ Re=16k, NCrit=9(ファイル {primary_name}):")
        print(f"    線形域揚力傾斜 dCL/dα = {slope:.4f} /deg")
        print(f"    設計点 CL_tail={CL_tail_d:+.3f} に対応する迎角 α ≈ {alpha_d:+.2f}°(線形範囲内)")
    else:
        print(f"  [警告] ポーラ未発見: {primary_path}")

    # ---------- 静安定余裕 ----------
    print("\n" + "─" * 72)
    print(" 静安定余裕 SM(縦の静安定の確認)")
    print("─" * 72)
    xnp_frac = neutral_point_frac(VH_TARGET, ar_h=AR_h)
    # SM は設計点(重心を翼弦 35%)・限界点(30%)の定義位置で評価:SM = x_NP/c − CG位置
    SM_d = xnp_frac - 0.35
    SM_l = xnp_frac - 0.30
    print(f"\n  中立点 x_NP/c = {xnp_frac:.3f}(AR_h={AR_h:.3f}, Vh={VH_TARGET})")
    print(f"  設計点(重心 翼弦35%): SM = x_NP/c − 0.35 = {SM_d*100:.1f} %")
    print(f"  限界点(重心 翼弦30%): SM = x_NP/c − 0.30 = {SM_l*100:.1f} %")

    # ---------- Step ③ 強度設計 ----------
    print("\n" + "─" * 72)
    print(" Step ③  強度設計(設計荷重 L_worst = 100 mN)")
    print("─" * 72)

    Ws_cfrp1 = section_modulus_circular(0.001)    # φ1 mm 翼内スパー
    Ws_cfrp2 = section_modulus_circular(0.002)    # φ2 mm ブーム

    # 水平尾翼
    M_h = root_moment_center_supported(L_WORST_DESIGN, B_BOOM)
    sigma_h = M_h / Ws_cfrp1
    SF_h = SIGMA_ULT_CFRP / sigma_h
    # EPS 単体(NACA0008 を等価矩形 0.5c × 0.08c で近似)
    W_eps_h = section_modulus_rectangle(0.5 * C_TAIL_H, 0.08 * C_TAIL_H)
    SF_eps_h = SIGMA_YIELD_EPS / (M_h / W_eps_h)
    print(f"\n  [水平尾翼]  M_root = L·b/8 = {M_h * 1000:.2f} mN·m")
    print(f"    φ1mm CFRP スパー: σ = {sigma_h / 1e6:.1f} MPa → SF = {SF_h:.0f}")
    print(f"    EPS 単体(補強なし): SF = {SF_eps_h:.1f}(二重安全性の裏付け)")

    # 垂直尾翼(U字 1枚あたり)
    M_v = root_moment_center_supported(L_WORST_DESIGN, B_VTAIL)
    sigma_v = M_v / Ws_cfrp1
    SF_v = SIGMA_ULT_CFRP / sigma_v
    print(f"\n  [垂直尾翼 U字 1枚あたり]  M_root = L·b_v/8 = {M_v * 1000:.2f} mN·m")
    print(f"    φ1mm CFRP スパー: σ = {sigma_v / 1e6:.2f} MPa → SF = {SF_v:.0f}")

    # ツインブーム
    M_boom = (L_WORST_DESIGN / 2.0) * ELL_BOOM
    sigma_boom = M_boom / Ws_cfrp2
    SF_boom = SIGMA_ULT_CFRP / sigma_boom
    print(f"\n  [ツインブーム φ2mm 1本]  M = (L/2)·ℓ_boom = {M_boom * 1000:.1f} mN·m")
    print(f"    σ = {sigma_boom / 1e6:.1f} MPa → SF = {SF_boom:.0f}")

    # ---------- サマリー ----------
    print("\n" + "=" * 72)
    print(" 確定値サマリー(計算書 §3 提供物)")
    print("=" * 72)
    print(f"  尾翼力      L_t = {L_t_d*1000:.1f} mN(設計点) / 設計荷重 {L_WORST_DESIGN*1000:.0f} mN")
    print(f"  水平尾翼    {B_BOOM*1000:.0f} × {C_TAIL_H*1000:.0f} mm  Sh = {Sh_adopt*1e4:.2f} cm²  AR_h = {AR_h:.3f}")
    print(f"  垂直尾翼    {B_VTAIL*1000:.0f} × {C_VTAIL*1000:.0f} mm × {N_VTAIL}枚  Sv = {Sv_adopt*1e4:.1f} cm²")
    print(f"  容積比      Vh = {Vh_real:.3f}  Vv = {Vv_real:.4f}")
    print(f"  静安定余裕  SM = {SM_d*100:.1f} %(設計点) / {SM_l*100:.1f} %(限界点)")
    print(f"  安全率      水平 SF={SF_h:.0f}  垂直 SF={SF_v:.0f}  ブーム SF={SF_boom:.0f}")
    print("\n[完了]")


if __name__ == '__main__':
    main()
