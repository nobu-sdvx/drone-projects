"""
aircraft_layout_v3.py
======================

機体配置図(簡易三面図)― v3:ツインブーム + U字型尾翼 構成

対象 : 概略計算書 §3/§4/§7 の図(機体配置図)
作成 : 2026-05-17
役割 : 現行 v3 構成(ツインブーム + U字尾翼 + ロッドフェアリング)の
       上面図・側面図を1枚で描く。旧 tail_layout(十字尾翼)の差し替え用。
出力 : plots/aircraft_layout_v3.png

注記 : 日本語フォント非依存のためラベルは英語表記。寸法は mm(機首=原点)。
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

# =============================================================================
# 機体寸法(mm、機首=0)― §2/§3/§4/§7 の確定値
# =============================================================================
POWERUP_LEN   = 220.0
BODY_W        = 18.2     # POWERUP/フェアリング 幅
BODY_H        = 14.6     # 側面 高さ
ROD_FRONT, ROD_REAR = 60.2, 205.2     # 露出ロッド(フェアリング被覆区間)
WING_LE, WING_TE = 128.0, 188.0       # 主翼前後縁
WING_SPAN     = 360.0
BOOM_Y        = 74.0     # ブーム左右位置(中心から)
BOOM_FRONT, BOOM_REAR = 158.0, 305.0  # ブーム前端・後端
HT_LE, HT_TE  = 269.0, 305.0          # 水平尾翼 前後縁(弦 36)
HT_SPAN       = 148.0
VT_CHORD, VT_HEIGHT = 30.0, 40.0      # 垂直尾翼 弦・高さ
PROP_X        = 230.0    # プロペラ面(モーター x=220 直後)
PROP_D        = 36.0
MOTOR_Y       = 19.925   # モーター軸 左右位置(軸間 39.85)
X_CG          = 149.0    # 設計重心
WING_SPAR_X   = 146.0    # 主翼桁(弦の30%)

BLUE   = "#4a78b5"
GREEN  = "#6aa84f"
YELLOW = "#e0c84a"
GRAY   = "#b0b0b0"
ORANGE = "#e08a3c"


def draw_top(ax):
    """上面図。"""
    # POWERUP 本体
    ax.add_patch(Rectangle((0, -BODY_W/2), POWERUP_LEN, BODY_W,
                           facecolor=GRAY, edgecolor="k", lw=1, label="POWERUP 4.0"))
    # ロッドフェアリング(露出ロッド被覆区間)
    ax.add_patch(Rectangle((ROD_FRONT, -BODY_W/2), ROD_REAR-ROD_FRONT, BODY_W,
                           facecolor="none", edgecolor=ORANGE, lw=1.2, ls="--"))
    # 主翼
    ax.add_patch(Rectangle((WING_LE, -WING_SPAN/2), WING_TE-WING_LE, WING_SPAN,
                           facecolor=BLUE, alpha=0.45, edgecolor="k", lw=1,
                           label="Main wing (AG36)"))
    # 主翼桁
    ax.plot([WING_SPAR_X, WING_SPAR_X], [-WING_SPAN/2, WING_SPAN/2],
            color="k", lw=1.0, ls=":", label="Wing spar (30%c)")
    # ツインブーム
    for s in (+1, -1):
        ax.plot([BOOM_FRONT, BOOM_REAR], [s*BOOM_Y, s*BOOM_Y],
                color="k", lw=2.2,
                label="Twin boom" if s == +1 else None)
    # 水平尾翼
    ax.add_patch(Rectangle((HT_LE, -HT_SPAN/2), HT_TE-HT_LE, HT_SPAN,
                           facecolor=GREEN, alpha=0.6, edgecolor="k", lw=1,
                           label="Horizontal tail"))
    # プロペラ(2 基)
    for s in (+1, -1):
        ax.add_patch(Circle((PROP_X, s*MOTOR_Y), PROP_D/2,
                            facecolor="none", edgecolor=ORANGE, lw=1.4, ls="-",
                            label="Propeller (D36)" if s == +1 else None))
    # 重心
    ax.plot(X_CG, 0, marker="o", color="red", ms=9, label="CG (149.0)")
    ax.plot(X_CG, 0, marker="+", color="white", ms=7, mew=1.5)

    # 寸法注記
    ax.annotate("", xy=(WING_LE, -198), xytext=(WING_TE, -198),
                arrowprops=dict(arrowstyle="<->", color="k"))
    ax.text((WING_LE+WING_TE)/2, -214, "chord 60", ha="center", fontsize=8)
    ax.annotate("", xy=(0, 250), xytext=(305, 250),
                arrowprops=dict(arrowstyle="<->", color="k"))
    ax.text(152, 258, "overall length 305", ha="center", fontsize=8)
    ax.text(BOOM_REAR+6, BOOM_Y, "boom pitch 148", fontsize=8, va="center")
    ax.text(WING_TE+10, 110, "span\n360", fontsize=8, ha="left")

    ax.set_title("Top view  (twin-boom + U-tail, v3)", fontsize=11, fontweight="bold")
    ax.set_xlabel("x from nose [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_xlim(-15, 345)
    ax.set_ylim(-235, 275)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7.5)


def draw_side(ax):
    """側面図。"""
    # POWERUP 本体
    ax.add_patch(Rectangle((0, -BODY_H/2), POWERUP_LEN, BODY_H,
                           facecolor=GRAY, edgecolor="k", lw=1))
    # ロッドフェアリング(かまぼこ:上半楕円のイメージを矩形で簡略)
    ax.add_patch(Rectangle((ROD_FRONT, -7.3), ROD_REAR-ROD_FRONT, 14.6,
                           facecolor=ORANGE, alpha=0.25, edgecolor=ORANGE, lw=1.2))
    ax.text((ROD_FRONT+ROD_REAR)/2, 10, "rod fairing", ha="center",
            fontsize=7, color=ORANGE)
    # 主翼(側面では薄い翼厚)
    ax.add_patch(Rectangle((WING_LE, 7), WING_TE-WING_LE, 6,
                           facecolor=BLUE, alpha=0.6, edgecolor="k", lw=1))
    ax.text((WING_LE+WING_TE)/2, 16, "main wing", ha="center", fontsize=7)
    # ツインブーム(側面では重なって1本に見える)
    ax.plot([BOOM_FRONT, BOOM_REAR], [0, 0], color="k", lw=2.5)
    ax.text((BOOM_FRONT+BOOM_REAR)/2, -6, "twin boom", ha="center", fontsize=7)
    # ブーム埋め込み区間(主翼後縁から30mm)
    ax.plot([BOOM_FRONT, WING_TE], [0, 0], color="red", lw=3.0, alpha=0.7)
    ax.text(195, -9, "boom embed 30 mm", fontsize=7, color="red", ha="left")
    # 水平尾翼
    ax.add_patch(Rectangle((HT_LE, 3), HT_TE-HT_LE, 3,
                           facecolor=GREEN, alpha=0.7, edgecolor="k", lw=1))
    # 垂直尾翼(U字:上向き)
    ax.add_patch(Rectangle((BOOM_REAR-VT_CHORD, 0), VT_CHORD, VT_HEIGHT,
                           facecolor=YELLOW, alpha=0.7, edgecolor="k", lw=1))
    ax.text(BOOM_REAR-VT_CHORD/2, VT_HEIGHT+4, "V-tail\n(upward)",
            ha="center", fontsize=7)
    # プロペラ(側面では円板を線で)
    ax.plot([PROP_X, PROP_X], [-PROP_D/2, PROP_D/2], color=ORANGE, lw=2.5)
    ax.text(PROP_X+4, PROP_D/2, "prop", fontsize=7, color=ORANGE)
    # 重心
    ax.plot(X_CG, 0, marker="o", color="red", ms=9)
    ax.plot(X_CG, 0, marker="+", color="white", ms=7, mew=1.5)
    ax.text(X_CG, -14, "CG", ha="center", fontsize=8, color="red")

    ax.set_title("Side view", fontsize=11, fontweight="bold")
    ax.set_xlabel("x from nose [mm]")
    ax.set_ylabel("z [mm]")
    ax.set_xlim(-15, 345)
    ax.set_ylim(-30, 60)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, "plots")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "aircraft_layout_v3.png")

    fig, (ax_top, ax_side) = plt.subplots(2, 1, figsize=(11, 8))
    draw_top(ax_top)
    draw_side(ax_side)
    fig.suptitle("Aircraft layout  v3  (W=36.1 g, twin-boom + U-tail)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  -> saved: {out_path}")


if __name__ == "__main__":
    main()
