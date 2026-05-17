# プロペラ設計 §5 — QPROP 環境 Mac→Windows 引き継ぎ書

**作成日**: 2026-05-14
**作成元**: Claude Code(macOS セッション)
**宛先**: Nobu(Windows 作業継続)
**状況**: Phase 1(環境構築 + 入力ファイル作成)を macOS で完了。Phase 2 以降を Windows で実施。

---

## 0. このドキュメントの目的

Windows で作業を再開する際、最初にこれを読めば「何が終わっていて、何を再現すればよいか」が分かるようにまとめる。

- §1 で **Mac で完了した作業内容**
- §2 で **リポジトリに残した成果物**
- §3 で **Windows でのセットアップ手順**
- §4 で **動作確認テスト**(Mac の結果と一致することを確認)
- §5 で **ハマりどころと回避策**
- §6 で **Phase 2 の次のステップ**(BEM Python 実装の入口)

---

## 1. Phase 1 で完了した作業(Mac、2026-05-14)

### 1.1 ツール選定の結論

| 候補 | 採否 | 理由 |
|---|---|---|
| **QPROP 1.22 + QMIL**(MIT Drela 公式) | **採用** | Fortran 77 で軽量、Mac/Windows どちらでもビルド可、BEM 解析 + MIL 設計の両方を一本で扱える |
| XROTOR(MIT 公式) | 不採用 | macOS で X11(XQuartz)依存が重い。QPROP+QMIL で機能カバー可能 |
| `xrotor-python`(PyPI、daniel-de-vries) | 不採用(Mac) | macOS arm64 + Python 3.14 で CMake ビルド失敗。**Windows + MinGW なら公式サポート環境**なので、Windows 移行後に再評価 |

### 1.2 macOS arm64 でビルドした成果

- **gfortran 15.2.0**(`brew install gcc` で導入)
- **QPROP 1.22 と QMIL** のバイナリ(`~/tools/Qprop/bin/qprop`、`~/tools/Qprop/bin/qmil`)
- 動作確認:README サンプル(CAM 6x3 + Speed-400)で正常動作

### 1.3 自設計プロペラの入力ファイル作成

設計書 v0.3 §3.4(Crazyflie 45mm を D=36mm にスケーリング)の数値を QPROP の `.prop` 形式に落とした。

### 1.4 自設計プロペラの設計点解析(Mac で実行済み)

```
$ ~/tools/Qprop/bin/qprop powerup_propeller.prop powerup_motor.mot 7.5 22100 0
```

| 量 | 値 | 設計書 v0.3 との対応 |
|---|---|---|
| 推力 T | **22.82 mN** | 必要 26.5 mN に対し 86%(§4.7「マージン薄い」と整合) |
| トルク Q | 0.157 mN·m | — |
| 軸出力 P_shaft | 0.363 W | — |
| プロペラ効率 η_prop | 0.472 | 改善対象 |
| 翼端マッハ M_tip | 0.118 | §4.5(40K rpm で 0.22)と整合 |
| 動作 Re(75%R 近傍) | 14,000-16,000 | §4.4(8k-20k 想定)と一致 |

**観察**:根本側(r/R < 0.30)で Cl が負 → β=29° と高い割に局所流入角が大きく、迎角がマイナス側に張り出している。設計改善の余地あり(根本捻り角を下げる)。

---

## 2. リポジトリに残した成果物

### 2.1 新規ファイル(本コミットで追加)

```
drone-projects/
├── qprop_data/
│   ├── powerup_propeller.prop      ← 自設計プロペラ(18 ステーション)
│   └── powerup_motor.mot           ← POWERUP モーター推定値
└── docs/
    └── propeller_qprop_setup_handoff_2026-05-14.md  ← 本ファイル
```

### 2.2 `powerup_propeller.prop` の中身概要

- **Nblades**: 2
- **R**: 0.018 m(= D/2 = 36mm/2)
- **翼型パラメータ**(NACA 4412 系の暫定値、**NeuralFoil 結果後に精緻化必須**):
  - CL0=0.50, CL_a=5.80
  - CLmin=-0.40, CLmax=1.10
  - CD0=0.025, CD2u=0.040, CD2l=0.025, CLCD0=0.50
  - REref=50000, REexp=-0.50
- **スケール係数**:Rfac=0.018(r/R 表記)、Cfac=0.001(mm→m)、Bfac=1.0(deg)
- **18 ステーション**(r/R = 0.15 → 1.00 を 0.05 刻み)

### 2.3 `powerup_motor.mot` の中身概要

- Type 1(brushed DC、簡易モデル)
- **Rmotor = 1.0 Ω**(マイクロモーター推定値)
- **Io = 0.05 A**(無負荷電流推定値)
- **Kv = 10,800 rpm/V**(40,000 rpm ÷ 3.7 V から逆算)

> ⚠️ モーターの実測パラメータは非公開。**RPM 指定で QPROP を呼ぶ場合、推力・トルク値はモーターパラメータに非依存**(電流・電圧の見積りだけが推定値ベース)。

---

## 3. Windows セットアップ手順

### 3.1 前提

- Windows 10/11
- Git for Windows 導入済(`git pull` が動く前提)
- Python 3.11 or 3.12 推奨(**3.14 は xrotor-python のビルドが通らないリスク高い**)

### 3.2 MSYS2 + MinGW で QPROP をビルド

| Step | 操作 |
|---|---|
| 1. MSYS2 導入 | https://www.msys2.org/ のインストーラをダウンロード・実行 |
| 2. ベース更新 | MSYS2 ターミナル起動 → `pacman -Syu` を2回ほど(再起動はさむ) |
| 3. ツールチェイン導入 | `pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran make` |
| 4. **MinGW 64-bit シェル** を起動(MSYS2 のショートカット内、青いアイコン) | これ以降のコマンドはこのシェルで実行 |
| 5. gfortran 確認 | `gfortran --version` で 13.x 以上が出れば OK |
| 6. QPROP 取得 | `cd ~ && mkdir -p tools && cd tools` <br>`curl -O https://web.mit.edu/drela/Public/web/qprop/qprop1.22.tar.gz` <br>`tar xzf qprop1.22.tar.gz` |
| 7. Makefile 編集 | `Qprop/bin/Makefile` を以下のように書き換え(Mac でやったのと同じ修正): <br>` FFLAGS = -O -fdefault-real-8 -std=legacy` <br>` FC = gfortran` <br>` FTNLIB =`(ifort 行はコメントアウト) |
| 8. ビルド | `cd Qprop/bin && make qprop qmil` |
| 9. 動作確認 | `cd ../runs && ../bin/qprop cam6x3 s400-6v-dd 0 0 8`(thrust=3.275 N 程度が出れば成功) |

**Windows では `LIBRARY_PATH` の SDK 回避策は不要**(あれは macOS 固有問題)。

### 3.3 自設計プロペラの動作確認

```
$ cd <repo>/qprop_data
$ ~/tools/Qprop/bin/qprop powerup_propeller.prop powerup_motor.mot 7.5 22100 0
```

→ **T = 22.82 mN, η_prop = 0.472**(Mac の結果)と **一致** すれば移行成功。

### 3.4 Python 環境

```
$ cd <repo>
$ python -m venv .venv
$ .venv\Scripts\activate    # PowerShell の場合
$ pip install numpy scipy matplotlib pandas neuralfoil
```

### 3.5 xrotor-python の導入(オプション、Phase 2 で必要なら)

```
$ pip install xrotor
```

- **公式サポート環境は Windows + MinGW**(macOS では失敗した)。
- 上記 §3.2 で MinGW gcc + gfortran が PATH 上にある状態で `pip install xrotor` を打つと、CMake → MinGW でビルドされる。
- 失敗したら `distutils.cfg` の compiler=mingw32 設定や、Python バージョン(3.11/3.12 推奨)を見直す。

---

## 4. 動作確認テスト

Windows セットアップ完了後、以下が **全て** 通れば Phase 1 が再現できている。

| # | テスト | 期待値 |
|---|---|---|
| 1 | `gfortran --version` | 13.x 以上 |
| 2 | `qprop` 単独実行 | usage メッセージ表示 |
| 3 | `qprop cam6x3 s400-6v-dd 0 0 8` | T ≈ 3.275 N, RPM ≈ 14,020 |
| 4 | `qprop powerup_propeller.prop powerup_motor.mot 7.5 22100 0` | T ≈ 22.82 mN, η ≈ 0.47 |
| 5 | `python -c "import numpy; print(numpy.__version__)"` | 2.x |

---

## 5. ハマりどころ(Mac で踏んだ罠の備忘)

| 症状 | 原因 | 対策(Windows ではどうか) |
|---|---|---|
| `make` で `ld: library 'System' not found` | macOS の SDK 探索パスが Homebrew gfortran から見えない | **Windows では発生しない**(MinGW が独自 libc を持っているため) |
| `pip install xrotor` で CMake エラー | macOS arm64 + Python 3.14 が未サポート | **Windows + MinGW + Python 3.11/3.12 なら公式環境** |
| Makefile デフォルトが `ifort` | Drela が Intel Fortran 想定 | **Windows でも同じ修正が必要**(`FC = gfortran` 化) |

---

## 6. Phase 2 の次のステップ

### 6.1 残作業の優先順位

1. **(優先1)NeuralFoil で翼型ポーラ生成** → `.prop` の翼型パラメータ(CL0, CL_a, CD0 等)を実値に更新
   - 候補:S1223、E205、薄キャンバ平板、AG12(設計書 v0.3 §4.4)
2. **(優先2)`propeller_design_analysis.py` 実装**(設計書 v0.3 §8.1)
   - Python の自作 BEM
   - QPROP を `subprocess` で呼んで検算比較する関数を含める
3. **(優先3)`propeller_structural_analysis.py` 実装**(設計書 v0.3 §8.2)
4. **(優先4)`uiuc_propeller_comparison.py` 実装**(設計書 v0.3 §8.3)
5. **(優先5)Phase 2.5**:Nobu 実測ブレード質量で `M_BLADE_TENTATIVE` を実値更新

### 6.2 Python ↔ QPROP 連携の実装ヒント

```python
import subprocess
from pathlib import Path

QPROP_BIN = Path.home() / "tools" / "Qprop" / "bin" / "qprop"  # Windows なら .exe

def run_qprop(prop_file: Path, motor_file: Path, V: float, rpm: float, volt: float = 0) -> dict:
    """
    QPROP を呼んで標準出力をパース。
    Returns: dict with 'thrust', 'torque', 'eff_prop', 'rpm_actual', etc.
    """
    result = subprocess.run(
        [str(QPROP_BIN), str(prop_file), str(motor_file), str(V), str(rpm), str(volt)],
        capture_output=True, text=True, check=True
    )
    # 出力行: "  V(m/s)    rpm    Dbeta   T(N)   ..." の下に1行データ
    # コメント行(#始まり)を飛ばして数値行をパース
    ...
```

### 6.3 期日

- **5/15(金)夕方**:Phase 2 完了(3 つの Python スクリプト + プロット 11 枚)
- **5/16(土)**:計算書 §5 原稿
- **5/18(月)**:第2回審査会

---

## 7. ブランチ運用ルール(再掲)

CLAUDE.md §4.1 と user-request §6 に従う:

- ✅ 既存ブランチでの新規ファイル追加・コミット OK
- ❌ **新規ブランチを Claude Code が勝手に切らない**
- ❌ 既存ファイル(主翼・尾翼スクリプト)の大幅変更禁止
- 「ブランチを切るべき?」は Nobu に確認してから

---

## 8. 連絡事項

- 本書に書いていない問題が出たら、設計書 v0.3(`claude/user-request/Propeller_design_spec_v0.3_2026-05-13.md`)と実装指示書(`claude/user-request/user-request_propeller.md`)を参照
- それでも解けない場合は管理スレッド(別 Claude 会話)に報告

---

**以上。Windows で `git pull` → §3 のセットアップを実行 → §4 で動作確認テスト → §6 の Phase 2 へ。**