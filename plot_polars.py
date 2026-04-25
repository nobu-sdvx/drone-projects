"""
plot_polars.py
--------------
役割:
  drone_projects/ 内の全翼型を NeuralFoil で解析し、
  ポーラー図3種を PNG で保存する:
    1. CL vs α(揚力曲線)
    2. CL vs CD(Lilienthal線図 / ドラッグポーラー)
    3. L/D vs α(効率曲線)
  保存先: スクリプトと同じフォルダ / plots/ 以下

  計算書への貼り付け用の「見栄えの良い図」を作るのが目的。
"""

import neuralfoil as nf
import numpy as np
import matplotlib.pyplot as plt
import os

# ---------- 設定 ----------
Re_CRUISE = 35000           # 巡航中心域
ALPHAS = np.linspace(-4, 16, 41)   # -4°〜16°を0.5°刻み
MODEL = "xlarge"            # バランス重視の精度

HERE = os.path.dirname(os.path.abspath(__file__))
DAT_DIR = os.path.join(HERE, "airfoils")
OUT_DIR = os.path.join(HERE, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- 翼型収集 ----------
candidates = sorted([f for f in os.listdir(DAT_DIR) if f.endswith(".dat")])
print(f"解析対象: {len(candidates)}翼型 @ Re={Re_CRUISE:,}")

# ---------- 一括解析 ----------
data = {}
for af in candidates:
    out = nf.get_aero_from_dat_file(
        filename=os.path.join(DAT_DIR, af),
        alpha=ALPHAS,
        Re=Re_CRUISE,
        model_size=MODEL,
    )
    CL = np.asarray(out["CL"])
    CD = np.asarray(out["CD"])
    data[af] = {
        "CL": CL,
        "CD": CD,
        "LD": CL / np.where(CD > 1e-6, CD, 1e-6),
    }

# ---------- 色分け: Top4を目立たせる ----------
HIGHLIGHT = {"ag36.dat", "ag35.dat", "sd7037.dat", "e387.dat"}
def style_for(af):
    name = af.replace(".dat", "")
    if af in HIGHLIGHT:
        return {"label": name, "linewidth": 2.2, "alpha": 1.0}
    else:
        return {"label": name, "linewidth": 1.0, "alpha": 0.5, "linestyle": "--"}

# ---------- Plot 1: CL vs α ----------
plt.figure(figsize=(9, 6))
for af in candidates:
    plt.plot(ALPHAS, data[af]["CL"], **style_for(af))
plt.axhline(0, color="gray", linewidth=0.5)
plt.axvline(0, color="gray", linewidth=0.5)
plt.xlabel("Angle of attack α [deg]", fontsize=11)
plt.ylabel("Lift coefficient $C_L$", fontsize=11)
plt.title(f"Lift curve (Re = {Re_CRUISE:,}, NeuralFoil {MODEL})", fontsize=12)
plt.legend(loc="lower right", fontsize=9, ncol=2)
plt.grid(True, alpha=0.3)
plt.tight_layout()
f1 = os.path.join(OUT_DIR, f"polar_CL_alpha_Re{Re_CRUISE}.png")
plt.savefig(f1, dpi=150)
plt.close()
print(f"  saved: {f1}")

# ---------- Plot 2: CL vs CD (ドラッグポーラー) ----------
plt.figure(figsize=(8, 7))
for af in candidates:
    plt.plot(data[af]["CD"], data[af]["CL"], **style_for(af))
plt.xlabel("Drag coefficient $C_D$", fontsize=11)
plt.ylabel("Lift coefficient $C_L$", fontsize=11)
plt.title(f"Drag polar (Re = {Re_CRUISE:,}, NeuralFoil {MODEL})", fontsize=12)
plt.legend(loc="lower right", fontsize=9, ncol=2)
plt.grid(True, alpha=0.3)
plt.xlim(0, 0.10)
plt.ylim(-0.5, 2.0)
plt.tight_layout()
f2 = os.path.join(OUT_DIR, f"polar_CL_CD_Re{Re_CRUISE}.png")
plt.savefig(f2, dpi=150)
plt.close()
print(f"  saved: {f2}")

# ---------- Plot 3: L/D vs α (効率曲線) ----------
plt.figure(figsize=(9, 6))
for af in candidates:
    plt.plot(ALPHAS, data[af]["LD"], **style_for(af))
plt.axhline(0, color="gray", linewidth=0.5)
plt.xlabel("Angle of attack α [deg]", fontsize=11)
plt.ylabel("Lift-to-drag ratio L/D", fontsize=11)
plt.title(f"Efficiency curve (Re = {Re_CRUISE:,}, NeuralFoil {MODEL})", fontsize=12)
plt.legend(loc="lower right", fontsize=9, ncol=2)
plt.grid(True, alpha=0.3)
plt.ylim(-5, 50)
plt.tight_layout()
f3 = os.path.join(OUT_DIR, f"polar_LD_alpha_Re{Re_CRUISE}.png")
plt.savefig(f3, dpi=150)
plt.close()
print(f"  saved: {f3}")

print(f"\n全3図を plots/ に保存しました。")
print(f"太い実線 = Top4候補 (AG36, AG35, SD7037, E387)")
print(f"細い破線 = 脱落組 (S1223, NACA4412, MH60, SD7003)")