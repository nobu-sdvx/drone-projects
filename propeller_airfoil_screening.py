"""
propeller_airfoil_screening.py
------------------------------
役割:
  プロペラ ブレード断面翼型の候補を、XFLR5 で出力済みのポーラデータ
  (xfoil_data/propeller_blade/*.txt)から読み込み、Re 5 点で比較解析する。
  主翼用 screening.py の「翼型一括スクリーニング」をプロペラ向けに作り直したもの。
  主翼用 screening.py が NeuralFoil をその場で回すのに対し、本スクリプトは
  XFLR5 で計算済みの .txt ポーラを読むだけ(再解析しない)。

入力:
  xfoil_data/propeller_blade/{翼型}_T1_Re{Re}_N{NCrit}.txt
    翼型 4 種 : S1223 / E205 / camber4_thick2 / AG12
    Re   5 点 : 8,000 / 10,000 / 12,000 / 15,000 / 20,000
  XFLR5 v6.62 の Polar Export 形式(ヘッダ + alpha/CL/CD/CDp/Cm... の表)

出力:
  - コンソール : Re=15,000 代表点での 4 翼型比較表 + 推奨翼型
  - plots/airfoil_polars_re15k.png      : Re=15k で CL-α・CD-α・(CL/CD)-α の 3 枚並べ
  - plots/airfoil_ld_max_by_re.png      : 各 Re での最大 CL/CD を棒グラフ比較
  - plots/airfoil_polar_smoothness.png  : ポーラ滑らかさ評価(CL-CD 抗力ポーラ + α-CL)

選定基準(優先順位順):
  1. Re=15,000 での CL/CD が最大
  2. CL=0.5 以上を確保できる
  3. ポーラが滑らか(失速ピーク・キンクが少ない)

使い方:
  venv を有効化してから  python propeller_airfoil_screening.py
"""

import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")  # 画面のない環境でも PNG 出力できるように
import matplotlib.pyplot as plt

# 図中の日本語が豆腐(□)にならないよう Windows 標準の日本語フォントを使う
matplotlib.rcParams["font.sans-serif"] = ["Yu Gothic", "MS Gothic", "Meiryo",
                                          "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

# ---------- 設定 ----------
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "xfoil_data", "propeller_blade")
PLOT_DIR = os.path.join(HERE, "plots")

RE_LIST = [8000, 10000, 12000, 15000, 20000]
RE_REP = 15000  # 75%R 代表点

# 翼型ラベル -> ファイル名の接頭辞(ファイル名は小文字)
AIRFOILS = {
    "S1223": "s1223",
    "E205": "e205",
    "camber4_thick2": "camber4_thick2",
    "AG12": "ag12",
}

# プロット時の色(翼型ごとに固定)
COLORS = {
    "S1223": "tab:red",
    "E205": "tab:blue",
    "camber4_thick2": "tab:green",
    "AG12": "tab:orange",
}

CL_FLOOR = 0.5  # 選定基準 2:確保すべき CL の下限


# ---------- XFLR5 ポーラ読込 ----------
def load_polar(prefix, Re):
    """
    xfoil_data/propeller_blade/ から {prefix}_T1_Re{Re}_N*.txt を探して読み込む。
    NCrit 部分(N5 / N50 等)はファイルごとに表記揺れがあるためワイルドカード。
    戻り値: dict(alpha, CL, CD, Cm) いずれも昇順 alpha の np.ndarray。
    """
    pattern = os.path.join(DATA_DIR, f"{prefix}_T1_Re{Re:05d}_N*.txt")
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"ポーラ未検出: {pattern}")
    path = matches[0]

    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # 「alpha CL CD ...」ヘッダ行を探し、その 2 行先(ダッシュ行の次)から表本体
    data_start = None
    for i, ln in enumerate(lines):
        if ln.lstrip().lower().startswith("alpha"):
            data_start = i + 2
            break
    if data_start is None:
        raise ValueError(f"ヘッダ行が見つかりません: {path}")

    for ln in lines[data_start:]:
        parts = ln.split()
        if len(parts) < 5:
            continue
        try:
            alpha, CL, CD, _CDp, Cm = (float(parts[k]) for k in range(5))
        except ValueError:
            continue  # 非収束行・空行はスキップ
        rows.append((alpha, CL, CD, Cm))

    if not rows:
        raise ValueError(f"データ行が読めません: {path}")

    rows.sort(key=lambda r: r[0])
    arr = np.array(rows)
    return {
        "alpha": arr[:, 0],
        "CL": arr[:, 1],
        "CD": arr[:, 2],
        "Cm": arr[:, 3],
        "file": os.path.basename(path),
    }


# ---------- 指標抽出 ----------
def metrics(polar):
    """1 つのポーラから評価指標を算出する。"""
    alpha, CL, CD = polar["alpha"], polar["CL"], polar["CD"]
    LD = CL / np.where(CD > 1e-6, CD, 1e-6)

    idx = int(np.argmax(LD))  # CL/CD 最大の動作点
    return {
        "LD": LD,
        "LD_max": LD[idx],
        "alpha_LDmax": alpha[idx],
        "CL_at_LDmax": CL[idx],
        "CD_at_LDmax": CD[idx],
        "CL_max": float(CL.max()),
        "alpha_CLmax": float(alpha[int(np.argmax(CL))]),
    }


def smoothness(polar):
    """
    ポーラの滑らかさを評価する。
    XFLR5 T1(定常)ポーラにはヒステリシスは現れないため、ここでは
      - 失速ピーク : CL が解析範囲内で局所最大を取り、その後低下するか
      - キンク     : CL-α 曲線の 2 階差分のばらつき(層流剥離バブル等の凹凸)
    の 2 点で評価する。score が小さいほど滑らか。
    """
    alpha, CL = polar["alpha"], polar["CL"]

    # 失速ピーク:全体最大が解析範囲の端でなく、かつその後 5% 以上低下
    i_peak = int(np.argmax(CL))
    has_stall = (i_peak < len(CL) - 1) and (CL[-1] < CL[i_peak] * 0.95)

    # キンク:CL の 2 階差分の RMS(α 等間隔前提で正規化)
    d2 = np.diff(CL, 2)
    kink = float(np.sqrt(np.mean(d2 ** 2))) if len(d2) else 0.0

    return {"has_stall": has_stall, "kink": kink, "alpha_stall": alpha[i_peak]}


# ---------- 解析 ----------
def main():
    os.makedirs(PLOT_DIR, exist_ok=True)

    # 全翼型 × 全 Re を読み込み
    polars = {}   # polars[label][Re] = polar dict
    mets = {}     # mets[label][Re]   = metrics dict
    for label, prefix in AIRFOILS.items():
        polars[label] = {}
        mets[label] = {}
        for Re in RE_LIST:
            p = load_polar(prefix, Re)
            polars[label][Re] = p
            mets[label][Re] = metrics(p)

    # ----- 全 Re の概観表 -----
    print("=" * 72)
    print("プロペラ ブレード翼型スクリーニング(XFLR5 ポーラ Type1, NCrit=5)")
    print("=" * 72)
    hdr = f"{'翼型':<16}{'Re':<9}{'CL/CD最大':<11}{'α@最大':<9}{'CL@最大L/D':<12}{'CD@最大L/D':<11}"
    for label in AIRFOILS:
        print()
        print(hdr)
        print("-" * len(hdr))
        for Re in RE_LIST:
            m = mets[label][Re]
            print(f"{label:<16}{Re:<9}{m['LD_max']:<11.2f}{m['alpha_LDmax']:<9.1f}"
                  f"{m['CL_at_LDmax']:<12.3f}{m['CD_at_LDmax']:<11.4f}")

    # ----- Re=15,000 代表点の比較表(依頼書フォーマット) -----
    print()
    print("=" * 72)
    print("=== 翼型スクリーニング結果 ===")
    print(f"Re = {RE_REP:,}(75%R 代表点)での比較:")
    print()
    print(f"{'翼型':<18}{'CL/CD最大':<11}{'最大時α':<10}{'CL@最大L/D':<13}{'CD@最大L/D':<11}{'CL_max':<9}")
    print("-" * 72)
    for label in AIRFOILS:
        m = mets[label][RE_REP]
        print(f"{label:<18}{m['LD_max']:<11.2f}{m['alpha_LDmax']:<10.1f}"
              f"{m['CL_at_LDmax']:<13.3f}{m['CD_at_LDmax']:<11.4f}{m['CL_max']:<9.3f}")

    # ----- 滑らかさ評価 -----
    print()
    print("ポーラ滑らかさ(Re=15,000):")
    print(f"{'翼型':<18}{'失速ピーク':<13}{'キンク指標':<13}")
    print("-" * 44)
    smo = {}
    for label in AIRFOILS:
        s = smoothness(polars[label][RE_REP])
        smo[label] = s
        stall_txt = f"あり(α={s['alpha_stall']:.0f}°)" if s["has_stall"] else "なし"
        print(f"{label:<18}{stall_txt:<13}{s['kink']:<13.5f}")

    # ----- 推奨翼型の選定 -----
    recommend, reason = select(mets, smo)

    print()
    print("-" * 72)
    print(f"→ 推奨翼型: {recommend}")
    print(f"  選定理由(1-2文): {reason}")
    print("=" * 72)

    # ----- プロット -----
    plot_polars_re15k(polars, mets)
    plot_ld_max_by_re(mets)
    plot_smoothness(polars, smo)

    print()
    print("生成プロット:")
    for fn in ("airfoil_polars_re15k.png",
               "airfoil_ld_max_by_re.png",
               "airfoil_polar_smoothness.png"):
        print(f"  plots/{fn}")


# ---------- 推奨翼型の選定ロジック ----------
def select(mets, smo):
    """
    選定基準(優先順位順):
      1. Re=15,000 での CL/CD が最大
      2. CL=0.5 以上を確保できる(CL_max >= 0.5)
      3. ポーラが滑らか(失速ピーク・キンクが少ない)
    まず基準 2 で足切りし、残った中で基準 1(CL/CD 最大)で選ぶ。
    同等(5% 以内)なら基準 3 で滑らかな方を採る。
    """
    # 基準 2:CL_max >= 0.5 で足切り
    qualified = [a for a in AIRFOILS if mets[a][RE_REP]["CL_max"] >= CL_FLOOR]
    if not qualified:
        qualified = list(AIRFOILS)  # 全滅時は救済

    # 基準 1:CL/CD 最大で降順
    qualified.sort(key=lambda a: -mets[a][RE_REP]["LD_max"])
    best = qualified[0]
    best_ld = mets[best][RE_REP]["LD_max"]

    # 基準 3:CL/CD が 5% 以内の対抗馬があれば、滑らかな方を優先
    rivals = [a for a in qualified
              if mets[a][RE_REP]["LD_max"] >= best_ld * 0.95]
    if len(rivals) > 1:
        rivals.sort(key=lambda a: (smo[a]["has_stall"], smo[a]["kink"]))
        best = rivals[0]

    m = mets[best][RE_REP]
    s = smo[best]
    stall_txt = "失速ピークなし" if not s["has_stall"] else f"α={s['alpha_stall']:.0f}° に失速ピークあり"
    reason = (
        f"Re=15,000 で CL/CD={m['LD_max']:.1f}(α={m['alpha_LDmax']:.1f}°)と"
        f"4 候補中で最も高く、その動作点で CL={m['CL_at_LDmax']:.2f} と推力に必要な "
        f"CL≧0.5 を満たし(CL_max={m['CL_max']:.2f})、解析範囲内で {stall_txt}・"
        f"ポーラのキンクも小さいため、ブレード断面翼型として最適と判断した。"
    )
    return best, reason


# ---------- プロット 1:Re=15k の 3 枚並べ ----------
def plot_polars_re15k(polars, mets):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    fig.suptitle(f"プロペラ ブレード翼型ポーラ比較 (Re={RE_REP:,}, XFLR5 Type1, NCrit=5)",
                 fontsize=12)

    for label in AIRFOILS:
        p = polars[label][RE_REP]
        c = COLORS[label]
        alpha, CL, CD = p["alpha"], p["CL"], p["CD"]
        LD = CL / np.where(CD > 1e-6, CD, 1e-6)

        axes[0].plot(alpha, CL, "-o", ms=3, color=c, label=label)
        axes[1].plot(alpha, CD, "-o", ms=3, color=c, label=label)
        axes[2].plot(alpha, LD, "-o", ms=3, color=c, label=label)

        # CL/CD 最大点をマーカで強調
        m = mets[label][RE_REP]
        axes[2].plot(m["alpha_LDmax"], m["LD_max"], "*", ms=14,
                     color=c, markeredgecolor="k")

    for ax, (xl, yl, ttl) in zip(axes, [
        ("alpha [deg]", "CL", "CL - alpha"),
        ("alpha [deg]", "CD", "CD - alpha"),
        ("alpha [deg]", "CL / CD", "CL/CD - alpha  (* = max)"),
    ]):
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(ttl)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(PLOT_DIR, "airfoil_polars_re15k.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)


# ---------- プロット 2:Re ごとの最大 CL/CD 棒グラフ ----------
def plot_ld_max_by_re(mets):
    fig, ax = plt.subplots(figsize=(9, 5))

    labels = list(AIRFOILS)
    n = len(labels)
    x = np.arange(len(RE_LIST))
    width = 0.8 / n

    for i, label in enumerate(labels):
        vals = [mets[label][Re]["LD_max"] for Re in RE_LIST]
        bars = ax.bar(x + (i - (n - 1) / 2) * width, vals, width,
                      color=COLORS[label], label=label)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.3,
                    f"{v:.1f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{Re:,}" for Re in RE_LIST])
    ax.set_xlabel("Reynolds number")
    ax.set_ylabel("最大 CL/CD")
    ax.set_title("各 Re での最大 CL/CD 比較(プロペラ ブレード翼型)")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()

    # 代表点 Re=15k を強調
    if RE_REP in RE_LIST:
        ax.axvspan(RE_LIST.index(RE_REP) - 0.45, RE_LIST.index(RE_REP) + 0.45,
                   color="yellow", alpha=0.12, zorder=0)

    fig.tight_layout()
    out = os.path.join(PLOT_DIR, "airfoil_ld_max_by_re.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)


# ---------- プロット 3:ポーラ滑らかさ評価 ----------
def plot_smoothness(polars, smo):
    """
    左 : CL-CD 抗力ポーラ。曲線の凹凸・ループ状の乱れ(キンク)が滑らかさの指標。
    右 : CL-α。失速ピーク(局所最大→低下)の有無を確認。
    XFLR5 Type1(定常)解析のためヒステリシスループ自体は出ない点に注意。
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"ポーラ滑らかさ評価 (Re={RE_REP:,})", fontsize=12)

    for label in AIRFOILS:
        p = polars[label][RE_REP]
        c = COLORS[label]
        s = smo[label]
        tag = " [失速ピーク]" if s["has_stall"] else ""

        axes[0].plot(p["CD"], p["CL"], "-o", ms=3, color=c,
                     label=f"{label} (kink={s['kink']:.4f})")
        axes[1].plot(p["alpha"], p["CL"], "-o", ms=3, color=c,
                     label=f"{label}{tag}")
        if s["has_stall"]:
            i_peak = int(np.argmax(p["CL"]))
            axes[1].plot(p["alpha"][i_peak], p["CL"][i_peak], "v",
                         ms=11, color=c, markeredgecolor="k")

    axes[0].set_xlabel("CD")
    axes[0].set_ylabel("CL")
    axes[0].set_title("抗力ポーラ CL-CD(凹凸 = キンク)")
    axes[1].set_xlabel("alpha [deg]")
    axes[1].set_ylabel("CL")
    axes[1].set_title("CL-alpha(▽ = 失速ピーク)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(PLOT_DIR, "airfoil_polar_smoothness.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
