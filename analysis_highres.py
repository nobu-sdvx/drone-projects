"""
analysis_highres.py
-------------------
役割:
  最終候補 Top3 (AG36, SD7037, E387) を NeuralFoil の最高精度モデル
  "xxxlarge" で再解析し、xlarge との差分を比較する。
  計算書に「感度解析(モデル精度依存性)を実施した」と書ける根拠になる。

出力:
  - コンソールに xlarge vs xxxlarge の数値比較表
  - plots/highres_comparison.png: CL-α を2モデル重ね描き
"""

import neuralfoil as nf
import numpy as np
import matplotlib.pyplot as plt
import os

# ---------- 設定 ----------
TOP3 = ["ag36.dat", "sd7037.dat", "e387.dat"]
Re_LIST = [24000, 32000, 50000]          # 巡航下限・巡航上限・最大速度域
ALPHAS = np.linspace(-2, 14, 33)

HERE = os.path.dirname(os.path.abspath(__file__))
DAT_DIR = os.path.join(HERE, "airfoils")
OUT_DIR = os.path.join(HERE, "plots")
os.makedirs(OUT_DIR, exist_ok=True)


def analyze(af, Re, model):
    out = nf.get_aero_from_dat_file(
        filename=os.path.join(DAT_DIR, af),
        alpha=ALPHAS,
        Re=Re,
        model_size=model,
    )
    CL = np.asarray(out["CL"])
    CD = np.asarray(out["CD"])
    LD = CL / np.where(CD > 1e-6, CD, 1e-6)
    idx = int(np.argmax(LD))
    return {
        "CL": CL, "CD": CD, "LD": LD,
        "LD_max": LD[idx],
        "alpha_LDmax": ALPHAS[idx],
        "CL_max": float(CL.max()),
        "CL_at_LDmax": CL[idx],
    }


# ---------- 比較表出力 ----------
print(f"{'Airfoil':<12}{'Re':<8}{'model':<10}{'L/D_max':<10}{'α_opt':<8}{'CL_max':<8}")
print("-" * 60)

results = {}
for af in TOP3:
    results[af] = {}
    for Re in Re_LIST:
        for model in ["xlarge", "xxxlarge"]:
            r = analyze(af, Re, model)
            results[af][(Re, model)] = r
            print(f"{af:<12}{Re:<8}{model:<10}"
                  f"{r['LD_max']:<10.2f}{r['alpha_LDmax']:<8.1f}{r['CL_max']:<8.2f}")
        # xlarge vs xxxlarge 差分
        rx  = results[af][(Re, "xlarge")]
        rxxx = results[af][(Re, "xxxlarge")]
        diff_ld = (rxxx["LD_max"] - rx["LD_max"]) / rx["LD_max"] * 100
        diff_cl = (rxxx["CL_max"] - rx["CL_max"]) / rx["CL_max"] * 100
        print(f"  → xxxlarge vs xlarge: ΔL/D={diff_ld:+.2f}%, ΔCL_max={diff_cl:+.2f}%")
    print()


# ---------- 可視化: xxxlarge での CL-α 曲線(Re=32000) ----------
Re_PLOT = 32000
plt.figure(figsize=(9, 6))
colors = {"ag36.dat": "tab:orange", "sd7037.dat": "tab:gray", "e387.dat": "tab:green"}
for af in TOP3:
    r_x  = results[af][(Re_PLOT, "xlarge")]
    r_xx = results[af][(Re_PLOT, "xxxlarge")]
    name = af.replace(".dat", "")
    c = colors[af]
    plt.plot(ALPHAS, r_x["CL"],  color=c, linestyle="--", linewidth=1.2,
             label=f"{name} (xlarge)", alpha=0.6)
    plt.plot(ALPHAS, r_xx["CL"], color=c, linestyle="-",  linewidth=2.2,
             label=f"{name} (xxxlarge)")
plt.axhline(0, color="gray", linewidth=0.5)
plt.axvline(0, color="gray", linewidth=0.5)
plt.xlabel("Angle of attack α [deg]", fontsize=11)
plt.ylabel("Lift coefficient $C_L$", fontsize=11)
plt.title(f"High-precision comparison: xlarge vs xxxlarge (Re = {Re_PLOT:,})",
          fontsize=12)
plt.legend(loc="lower right", fontsize=9)
plt.grid(True, alpha=0.3)
plt.tight_layout()
f = os.path.join(OUT_DIR, "highres_comparison.png")
plt.savefig(f, dpi=150)
plt.close()
print(f"保存: {f}")
print("\n読み解き: 破線(xlarge)と実線(xxxlarge)がほぼ重なっていれば")
print("        xlargeモデルの精度で十分ということ。差が大きい=注意。")