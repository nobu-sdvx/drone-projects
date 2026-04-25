"""
flight_envelope.py
------------------
役割:
  飛行包絡線を可視化する。
  「ある速度Vで飛ぶために必要なCL」と
  「その速度での各翼型の出せる最大CL」を同じ図に描き、
  「どの速度域で飛べるか/飛べないか」を定量評価する。

原理:
  水平定常飛行では L = W なので、
    CL_required = 2·W / (ρ·V²·S)
  これが翼型のCL_max以下なら飛べる。

出力:
  - plots/flight_envelope.png: 要求CLと利用可能CLの比較
  - コンソール: 各翼型の失速速度(最低飛行速度)
"""

import neuralfoil as nf
import numpy as np
import matplotlib.pyplot as plt
import os

# ---------- 機体諸元(要件定義書より) ----------
W_target   = 0.325 * 9.81 / 1000 * 1000   # 見込み 32.5g → N  ... 単位を素直に
W_target   = 0.0325 * 9.81    # 32.5 g → 0.319 N
W_max      = 0.040  * 9.81    # 40 g(上限) → 0.392 N
S_wing     = 0.0216           # 翼面積 216 cm² → 0.0216 m²
CHORD      = 0.060            # 翼弦 60 mm → 0.060 m
RHO        = 1.225            # 空気密度 [kg/m³]
NU         = 1.5e-5           # 空気動粘性係数 [m²/s]

# ---------- 評価対象 ----------
TOP3 = ["ag36.dat", "sd7037.dat", "e387.dat"]
COLORS = {"ag36.dat": "tab:orange", "sd7037.dat": "tab:gray", "e387.dat": "tab:green"}

# 飛行速度スイープ
V_sweep = np.linspace(3.0, 12.0, 37)   # 3〜12 m/s

HERE = os.path.dirname(os.path.abspath(__file__))
DAT_DIR = os.path.join(HERE, "airfoils")
OUT_DIR = os.path.join(HERE, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- 必要CL計算 ----------
def CL_required(V, W):
    return 2.0 * W / (RHO * V**2 * S_wing)

CL_req_target = CL_required(V_sweep, W_target)   # 32.5g
CL_req_max    = CL_required(V_sweep, W_max)      # 40g

# ---------- 各翼型の CL_max(Re) を取得 ----------
# 速度V毎に Re を計算して NeuralFoil で CL_max を求める
def Re_from_V(V):
    return V * CHORD / NU

alphas = np.linspace(-2, 16, 37)
CL_max_available = {af: [] for af in TOP3}

print("CL_max 計算中...")
for V in V_sweep:
    Re = Re_from_V(V)
    for af in TOP3:
        out = nf.get_aero_from_dat_file(
            filename=os.path.join(DAT_DIR, af),
            alpha=alphas,
            Re=Re,
            model_size="xxxlarge",      # せっかくなので最高精度
        )
        CL_max_available[af].append(float(np.asarray(out["CL"]).max()))

for af in TOP3:
    CL_max_available[af] = np.array(CL_max_available[af])

# ---------- 失速速度の特定(CL_required == CL_max_available の交点) ----------
print("\n失速速度(= 最低飛行速度):")
print(f"{'Airfoil':<10}{'V_stall(32.5g)':<18}{'V_stall(40g)':<15}")
print("-" * 45)
V_stall_records = {}
for af in TOP3:
    # 32.5g の場合
    # V_stall は CL_required(V) == CL_max(V) となる V
    diff_t = CL_max_available[af] - CL_req_target
    diff_m = CL_max_available[af] - CL_req_max

    def find_crossover(diff, V):
        # 低速側から見て、最初に diff が正になる場所を失速速度とする
        # (速度が高ければ必要CLが小さいので diff は正、速度が低いと負)
        for i in range(1, len(V)):
            if diff[i - 1] < 0 <= diff[i]:
                # 線形補間
                frac = -diff[i - 1] / (diff[i] - diff[i - 1])
                return V[i - 1] + frac * (V[i] - V[i - 1])
        return None

    vs_t = find_crossover(diff_t, V_sweep)
    vs_m = find_crossover(diff_m, V_sweep)
    V_stall_records[af] = (vs_t, vs_m)
    name = af.replace(".dat", "").upper()
    vs_t_str = f"{vs_t:.2f} m/s" if vs_t else "飛行不可"
    vs_m_str = f"{vs_m:.2f} m/s" if vs_m else "飛行不可"
    print(f"{name:<10}{vs_t_str:<18}{vs_m_str:<15}")

# ---------- プロット ----------
fig, ax = plt.subplots(figsize=(11, 7))

# 要求CL曲線(負の領域=現実的に困難なライン)
ax.plot(V_sweep, CL_req_target, "k-",  linewidth=2.0, label="Required $C_L$ (32.5g)")
ax.plot(V_sweep, CL_req_max,    "k--", linewidth=1.5, label="Required $C_L$ (40g)")

# 各翼型の CL_max 曲線
for af in TOP3:
    name = af.replace(".dat", "").upper()
    ax.plot(V_sweep, CL_max_available[af], color=COLORS[af],
            linewidth=2.2, label=f"{name} CL_max")

# 巡航速度域の塗りつぶし
ax.axvspan(6, 8, alpha=0.15, color="green", label="Cruise 6–8 m/s")
ax.axvline(9, color="red", linestyle=":", linewidth=1.2, alpha=0.7, label="Max 9 m/s (PowerUp 4.0)")
ax.axvline(4, color="purple", linestyle=":", linewidth=1.2, alpha=0.7, label="Target stall 4 m/s")

ax.set_xlabel("Flight speed V [m/s]", fontsize=11)
ax.set_ylabel("Lift coefficient $C_L$", fontsize=11)
ax.set_title("Flight envelope: required vs available CL", fontsize=12)
ax.legend(loc="upper right", fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(3, 12)
ax.set_ylim(0, 2.0)

# 注釈
ax.annotate(
    "Flyable\n(Required < Available)",
    xy=(8, 0.5), fontsize=10, ha="center",
    bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.7)
)
ax.annotate(
    "Stall zone",
    xy=(4, 1.7), fontsize=10, ha="center", color="red",
    bbox=dict(boxstyle="round", facecolor="mistyrose", alpha=0.7)
)

plt.tight_layout()
f = os.path.join(OUT_DIR, "flight_envelope.png")
plt.savefig(f, dpi=150)
plt.close()
print(f"\n保存: {f}")

# ---------- 結論 ----------
print("\n" + "=" * 60)
print("結論")
print("=" * 60)
print(f"現在の設計(翼面積 {S_wing*1e4:.0f} cm²)では、")
print(f"要件『失速速度 4 m/s以下』の達成は物理的に困難。")
print(f"改善策の候補:")
print(f"  1. 翼面積を拡大: 例 216→300 cm² (翼幅 360→500mm)")
print(f"  2. 失速速度の目標を緩和: 例 4 → 5 m/s以下")
print(f"  3. 総重量を更に軽減: 40→32.5g厳守")
print(f"この判断は第2回審査会(5/18)の概略計算書までに要整合。")