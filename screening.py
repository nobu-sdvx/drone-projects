"""
screening.py
------------
役割:
  drone_projects/ フォルダ内の全ての.datファイル(翼型)を
  目標レイノルズ数域(Re=20k, 35k, 50k)で一気に解析し、
  L/D最大値とCL最大値を表にして出力する。

使い方:
  このファイルと同じ階層の airfoils/ フォルダに .dat ファイルを置いてから、
  python screening.py
"""

import neuralfoil as nf
import numpy as np
import os

# ---------- 設定 ----------
Re_list = [20000, 35000, 50000]              # 要件の飛行速度域
alphas = np.linspace(-2, 14, 33)             # -2°〜14°を0.5°刻み
MODEL = "xlarge"                              # xxsmall〜xxxlargeの8段階

# ---------- 候補自動収集 ----------
here = os.path.dirname(os.path.abspath(__file__))
DAT_DIR = os.path.join(here, "airfoils")
candidates = sorted([f for f in os.listdir(DAT_DIR) if f.endswith(".dat")])
print(f"検出した翼型: {len(candidates)}個")
print(f"→ {', '.join(candidates)}\n")

# ---------- 解析 ----------
header = f"{'Airfoil':<14}{'Re':<8}{'L/D_max':<10}{'α@LDmax':<10}{'CL_max':<10}{'CL@LDmax':<10}"
print(header)
print("-" * len(header))

results = {}  # 後で使うため保存

for af in candidates:
    results[af] = {}
    for Re in Re_list:
        out = nf.get_aero_from_dat_file(
            filename=os.path.join(DAT_DIR, af),
            alpha=alphas,
            Re=Re,
            model_size=MODEL,
        )
        CL, CD = out["CL"], out["CD"]
        LD = CL / np.where(CD > 1e-6, CD, 1e-6)  # ゼロ割防止
        idx = np.argmax(LD)

        results[af][Re] = {
            "LD_max": LD[idx],
            "alpha_LDmax": alphas[idx],
            "CL_max": CL.max(),
            "CL_at_LDmax": CL[idx],
        }

        print(f"{af:<14}{Re:<8}{LD[idx]:<10.1f}{alphas[idx]:<10.1f}"
              f"{CL.max():<10.2f}{CL[idx]:<10.2f}")
    print()

# ---------- 簡易ランキング(Re=35000基準) ----------
print("=" * 60)
print("Re=35,000 での L/D_max ランキング")
print("=" * 60)
ranked = sorted(candidates, key=lambda a: -results[a][35000]["LD_max"])
for i, af in enumerate(ranked, 1):
    r = results[af][35000]
    print(f"{i}. {af:<14} L/D={r['LD_max']:.1f}  CL_max={r['CL_max']:.2f}")