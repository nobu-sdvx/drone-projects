# Claude Code 向け §5 プロペラ設計 実装指示書

**作成日**: 2026-05-13(水)
**作成元**: §5 管理スレッド(Claude Sonnet)
**宛先**: Claude Code(別会話、コマンドライン経由)
**目的**: 設計書 v0.3 の数値根拠を Python スクリプトで実装し、計算書 §5 の数値・図表を生成
**期限**: 2026-05-15(金)夕方までに完成 → 5/16 計算書原稿 → 5/18 審査会

---

## 0. このドキュメントの使い方

このファイルだけ読めば、Claude Code が `drone-projects/` リポジトリ内で §5 プロペラ設計の Python 実装に着手できるようにまとめてあります。

- §1 で「開発環境とリポジトリ構成」を把握
- §2 で「凍結された設計値」を確認(これらは変更してはいけない)
- §3 で「実装すべき 3 つのスクリプト」を理解
- §4 で「Crazyflie UIUC データの活用方法」を確認
- §5 で「出力すべき数値・図表」を把握
- §6 で「ブランチ運用ルール」を厳守(**勝手にブランチを切らない**)
- §7 で「コーディング規約」を確認
- §8 で「自己レビュー観点」を確認のうえ、管理スレッドに戻る

---

## 1. 開発環境とリポジトリ構成

### 1.1 環境

- **OS**: macOS(Nobu の Mac、ホームディレクトリ `/Users/naitonobuhito`)
- **Python**: PyCharm + venv(`drone-projects/` 配下に venv 設定済み)
- **主要パッケージ**: numpy, scipy, matplotlib, pandas, neuralfoil
- **numpy 注意**: numpy 2.x で `np.trapz` 削除済み → **`np.trapezoid` を使用**(`master_todo_2026-05-07.md` で確認)

### 1.2 既存リポジトリ構成

```
drone-projects/
├── CLAUDE.md                          # Claude Code 向けブリーフィング(本ファイルと併用)
├── airfoils/                          # 翼型 .dat ファイル(主翼解析用)
├── plots/                             # 自動生成される図
├── xfoil_data/                        # XFLR5 .txt × 16(主翼解析)
│   └── propeller_blade/               # ★ §5 で新規追加(ブレード断面翼型用)
├── xflr5_projects/                    # .xfl ファイル
├── uiuc_propeller_data/               # ★ §5 で新規追加
│   └── volume-2/data/                 # Crazyflie 等の UIUC データ(Nobu が事前取得済み)
├── screening.py                       # NeuralFoil 翼型スクリーニング(主翼で使用、流用可)
├── analysis_highres.py                # XFLR5 ポーラ取得(流用可)
├── mainwing_structural_analysis.py    # 主翼構造解析 170行(構造解析の流儀の参考)
├── weight_cg_analysis.py              # 重量配分・重心解析
└── flight_envelope.py                 # 飛行包絡線
```

### 1.3 §5 で新規作成するファイル

```
drone-projects/
├── propeller_design_analysis.py       # ★ メイン: BEM、推力・効率
├── propeller_structural_analysis.py   # ★ 強度: 応力、SF、感度解析
├── uiuc_propeller_comparison.py       # ★ UIUC 比較: Crazyflie 等との照合
├── xfoil_data/propeller_blade/        # ★ ブレード断面翼型ポーラ(NeuralFoil 出力)
└── (Google Drive)
    └── ドローン/data関連/グラフ/
        └── 2026-05-15_propeller_design/   # ★ 生成図 PNG 保存先
```

### 1.4 ファイル保存ルール(主翼 §2 と同様)

- **ソースコード・データ**: GitHub の `drone-projects/` に保存
- **生成された図(PNG)**: 一旦 `plots/` に保存 → 完成後に手動で Google Drive にコピー(Nobu が実施)

---

## 2. 凍結された設計値(変更不可)

以下の数値は設計書 v0.3 で確定済み。**スクリプトでは入力定数として扱い、コード内でハードコーディングせず `dataclass` や `config` で管理** すること。

### 2.1 機体側の入力値

```python
# 機体物性
W_TOTAL = 0.348              # 機体総重量 [N](35.5 g)
V_CRUISE = 7.5               # 巡航速度 [m/s]
V_MAX = 9.0                  # 最大速度 [m/s]
RHO = 1.225                  # 空気密度 [kg/m³]
MU = 1.8e-5                  # 動粘度 [kg/m·s]
NU = 1.5e-5                  # 動粘度比(運動学的) [m²/s]
Q_CRUISE = 34.45             # 動圧 @ V=7.5 m/s [Pa]

# 主翼・尾翼から渡された値(凍結)
L_WING_DESIGN = 0.389        # 主翼揚力(設計点、W+|L_t|) [N]
L_WING_LIMIT = 0.398         # 主翼揚力(限界点) [N]
D_TOTAL_AIRFRAME = 0.0429    # 全機抗力(主翼+尾翼) [N]
D_FUSELAGE_EST = 0.010       # 胴体抗力推定(保守側) [N]
D_TOTAL = 0.0529             # 全機抗力(推定) [N]
T_REQUIRED_PER_MOTOR = D_TOTAL / 2   # 1モーター必要推力 = 0.0265 [N]
```

### 2.2 POWERUP 4.0 仕様

```python
# POWERUP モーター
N_MOTORS = 2
RPM_MAX = 40000              # 公式仕様 [rpm]
T_MAX_PER_MOTOR = 0.186      # 1モーター最大推力 [N](T/W=2 から)
AXIS_DISTANCE = 0.03985      # モーター軸間距離 [m](Nobu 実測)
SHAFT_DIAMETER = 0.001       # モーター軸径 [m](実測確認)
```

### 2.3 自作プロペラ設計値(凍結)

```python
# 自作プロペラ寸法
D_PROPELLER = 0.036          # 直径 [m]
P_PROPELLER = 0.027          # ピッチ [m]
P_OVER_D = P_PROPELLER / D_PROPELLER   # 0.75
N_BLADES = 2

# 動作点
RPM_CRUISE = 22100           # 巡航動作RPM(逆算) [rpm]
RPM_MAX_OP = 26600           # 最大動作RPM(V=9m/s) [rpm]
SLIP = 0.25                  # スリップ率

# ハブ仕様
HUB_HOLE_DIAMETER = 0.00095  # ハブ穴径 [m](軸径 -0.05mm の締めしろ)
HUB_LENGTH = 0.008           # ハブ長 [m]
HUB_OUTER_DIAMETER = 0.005   # ハブ外径 [m]

# ブレード根本断面(強度設計から逆算)
C_ROOT = 0.005               # 根本弦長 [m]
T_ROOT = 0.0015              # 根本厚さ [m]
```

### 2.4 ブレード形状(Crazyflie 45mm からスケーリング)

```python
# r/R, c/R, β [°]
BLADE_GEOMETRY = [
    (0.15, 0.245, 29.2),    # 根本
    (0.20, 0.263, 26.8),
    (0.25, 0.281, 24.6),
    (0.30, 0.297, 22.9),
    (0.35, 0.316, 21.5),
    (0.40, 0.333, 20.2),
    (0.45, 0.349, 19.0),
    (0.50, 0.363, 18.1),
    (0.55, 0.372, 17.5),
    (0.60, 0.379, 17.1),
    (0.65, 0.381, 16.8),
    (0.70, 0.382, 16.5),
    (0.75, 0.379, 16.1),    # 代表点
    (0.80, 0.372, 15.7),
    (0.85, 0.351, 15.6),
    (0.90, 0.311, 15.9),
    (0.95, 0.228, 15.5),
    (1.00, 0.046, 11.5),    # 先端
]
# 出典: cfnq_45p1_geom.txt(UIUC Volume 2)
# 我々の D=36mm にスケーリング: c [mm] = (c/R) × 18mm
```

### 2.5 材料定数(ABS)

```python
RHO_ABS = 1.04 * 1000        # 密度 [kg/m³]
SIGMA_Y_ABS = 40e6           # 降伏応力 [Pa]
E_ABS = 2.0e9                # 弾性率 [Pa]
TAU_EPOXY = 10e6             # エポキシせん断強度 [Pa](保守値)
ETA_EFFECTIVE = 0.10         # 接着有効率(保守値)
```

### 2.6 仮置き値(明日 5/14 確定予定)

```python
M_BLADE_TENTATIVE = 0.3e-3   # ブレード質量 [kg](0.3g 仮値)
# Nobu が 0.01g 精度の天秤で測定後、実測値に更新する
# 感度解析で 0.2-0.5g の範囲を評価
```

---

## 3. 実装すべき 3 つのスクリプト

### 3.1 `propeller_design_analysis.py`(メイン解析)

**目的**: BEM(ブレード要素運動量理論)で推力・効率曲線を計算し、UIUC データと比較する。

**期待する出力**:
1. プロペラ三面図(D・P・ハブ寸法、ブレード形状)
2. ブレード弦長分布 c(R) プロット(Crazyflie との重ね合わせ)
3. ブレードねじり角分布 β(R) プロット(Crazyflie との重ね合わせ)
4. 推力 vs RPM 曲線(V=0、5、7.5、9 m/s の 4 ケース)
5. 効率 η vs 前進比 J 曲線
6. UIUC Crazyflie 静推力との比較プロット

**主要な関数構成**:

```python
from dataclasses import dataclass
import numpy as np
from typing import Optional

@dataclass
class BladeGeometry:
    """プロペラブレードの形状を保持するデータクラス"""
    diameter: float          # 直径 [m]
    pitch: float             # ピッチ [m]
    n_blades: int            # ブレード数
    radial_stations: np.ndarray   # r/R 配列
    chord_distribution: np.ndarray  # c [m] 配列
    twist_distribution: np.ndarray  # β [rad] 配列
    airfoil_name: str        # 翼型名

def load_blade_geometry_from_crazyflie(D: float = 0.036) -> BladeGeometry:
    """
    Crazyflie 45mm の geom データを D=36mm にスケーリング。
    Returns: BladeGeometry インスタンス
    """
    pass

def compute_reynolds_distribution(geom: BladeGeometry, rpm: float, V: float) -> np.ndarray:
    """
    各半径での動作 Re を計算。
    Re = V_rel × c × ρ / μ
    V_rel = sqrt(V² + (ω·r)²)
    """
    pass

def load_airfoil_polars(airfoil_name: str, re_values: list) -> dict:
    """
    NeuralFoil で翼型ポーラ取得。
    各 Re ごとに CL-α、CD-α 曲線を返す。
    """
    pass

def bem_analysis(
    geom: BladeGeometry,
    rpm: float,
    V_inf: float,
    polars: dict,
    rho: float = 1.225
) -> dict:
    """
    BEM(ブレード要素運動量理論)で推力・トルク・効率を計算。
    
    BEM の基本式:
    dT = ρ·c·dr·(V_a² + V_t²)·(C_L·cosφ - C_D·sinφ)/2
    dQ = ρ·c·dr·r·(V_a² + V_t²)·(C_L·sinφ + C_D·cosφ)/2
    φ = atan(V_a / V_t)
    α = β - φ
    
    各半径で誘導速度を反復解法で収束させる。
    
    Returns:
        dict with keys: 'thrust', 'torque', 'efficiency', 'power', 'sectional_results'
    """
    pass

def compute_thrust_curve(geom: BladeGeometry, V_range: np.ndarray, rpm_range: np.ndarray) -> dict:
    """
    推力 vs V × RPM の 2D マップ生成。
    Returns: dict with 'V', 'RPM', 'T'(2次元配列)
    """
    pass

def compute_efficiency_curve(geom: BladeGeometry, rpm: float, V_range: np.ndarray) -> dict:
    """
    効率 η vs 前進比 J = V / (n·D) 曲線生成。
    """
    pass

def plot_geometry_3view(geom: BladeGeometry, save_path: str) -> None:
    """
    プロペラ三面図プロット。
    """
    pass

def plot_blade_distribution_comparison(self_geom: BladeGeometry, crazyflie_geom: BladeGeometry, save_path: str) -> None:
    """
    自作 vs Crazyflie のブレード c(R), β(R) 比較プロット。
    """
    pass

if __name__ == "__main__":
    # 1. ブレード形状の確定
    geom = load_blade_geometry_from_crazyflie(D=0.036)
    
    # 2. 動作Re域の確認
    rpm_design = 22100
    re_dist = compute_reynolds_distribution(geom, rpm_design, V_CRUISE)
    print(f"動作Re域: {re_dist.min():.0f} - {re_dist.max():.0f}")
    
    # 3. 翼型ポーラ取得
    polars = load_airfoil_polars(geom.airfoil_name, re_values=[8000, 12000, 16000, 20000])
    
    # 4. BEM 計算
    result = bem_analysis(geom, rpm_design, V_CRUISE, polars)
    print(f"巡航推力: {result['thrust']*1000:.1f} mN (目標 26.5 mN)")
    print(f"巡航効率: {result['efficiency']:.3f}")
    
    # 5. 推力曲線生成
    V_range = np.linspace(0, 10, 11)
    rpm_range = np.array([15000, 20000, 25000, 30000, 35000, 40000])
    thrust_map = compute_thrust_curve(geom, V_range, rpm_range)
    
    # 6. 効率曲線生成
    eff_curve = compute_efficiency_curve(geom, rpm_design, V_range)
    
    # 7. プロット生成
    plot_geometry_3view(geom, "plots/propeller_3view.png")
    # ... (他のプロット)
```

### 3.2 `propeller_structural_analysis.py`(強度解析)

**目的**: 設計したプロペラ形状に対し、応力・安全率を計算する。

**期待する出力**:
1. ブレード根応力分布プロット
2. SF 感度ヒートマップ(m_blade × RPM)
3. 軸-ハブ嵌合保持力 SF プロット

**主要な関数構成**:

```python
def compute_centrifugal_force(blade_mass: float, r_cg: float, omega: float) -> float:
    """
    遠心力 F = m × r × ω²
    
    Returns:
        F [N]
    """
    pass

def compute_blade_root_tensile_stress(F_centrifugal: float, area_root: float) -> float:
    """ σ_T = F / A """
    pass

def compute_blade_root_bending_stress(T_max: float, L_blade: float, c_root: float, t_root: float) -> float:
    """
    曲げ応力 σ_B = M / Z
    M = ½ × T × L(等分布荷重近似)
    Z = b × h² / 6
    """
    pass

def compute_safety_factor(sigma_combined: float, sigma_y: float = SIGMA_Y_ABS) -> float:
    """ SF = σ_y / σ_combined """
    pass

def compute_hub_holding_force(
    adhesive_area: float,
    tau: float = TAU_EPOXY,
    eta: float = ETA_EFFECTIVE
) -> float:
    """
    軸-ハブ嵌合の保持力 F_hold = τ × A × η
    A = π × shaft_d × hub_length
    """
    pass

def sensitivity_analysis_rpm_mass(
    geom: BladeGeometry,
    rpm_range: np.ndarray,
    mass_range: np.ndarray
) -> np.ndarray:
    """
    RPM × m_blade の 2D 感度解析。
    Returns: SF の 2次元配列
    """
    pass

def compute_tip_mach(omega: float, R: float, c_sound: float = 343) -> float:
    """ M = ω·R / c """
    pass

def plot_stress_distribution_radial(geom: BladeGeometry, rpm: float, save_path: str) -> None:
    """
    各半径での応力分布プロット。
    """
    pass

def plot_sf_sensitivity_heatmap(rpm_range, mass_range, sf_matrix, save_path: str) -> None:
    """
    SF 感度ヒートマップ(横軸 RPM、縦軸 m_blade、カラー SF)。
    SF=4 のコンター線を強調。
    """
    pass

if __name__ == "__main__":
    geom = load_blade_geometry_from_crazyflie(D=0.036)
    
    # 1. 三点評価(40K, 27K, 22K rpm)
    for rpm_name, rpm in [("40K最悪", 40000), ("27K最大動作", 26600), ("22K巡航", 22100)]:
        omega = 2 * np.pi * rpm / 60
        F_cf = compute_centrifugal_force(M_BLADE_TENTATIVE, 0.0135, omega)
        A_root = C_ROOT * T_ROOT
        sigma_T = compute_blade_root_tensile_stress(F_cf, A_root)
        sigma_B = compute_blade_root_bending_stress(T_MAX_PER_MOTOR, 0.016, C_ROOT, T_ROOT)
        sigma_total = sigma_T + sigma_B
        SF = compute_safety_factor(sigma_total)
        print(f"{rpm_name}({rpm}rpm): σ={sigma_total/1e6:.2f} MPa, SF={SF:.1f}")
    
    # 2. 軸-ハブ嵌合
    A_adh = np.pi * SHAFT_DIAMETER * HUB_LENGTH
    F_hold = compute_hub_holding_force(A_adh)
    SF_hub = F_hold / T_MAX_PER_MOTOR
    print(f"軸-ハブ嵌合 SF: {SF_hub:.0f}")
    
    # 3. 翼端マッハ
    omega_max = 2 * np.pi * 40000 / 60
    M_tip = compute_tip_mach(omega_max, D_PROPELLER / 2)
    print(f"翼端マッハ数(40K rpm): {M_tip:.3f}")
    
    # 4. 感度解析
    rpm_range = np.linspace(15000, 40000, 11)
    mass_range = np.linspace(0.0002, 0.0005, 7)
    sf_matrix = sensitivity_analysis_rpm_mass(geom, rpm_range, mass_range)
    plot_sf_sensitivity_heatmap(rpm_range, mass_range, sf_matrix, "plots/sf_heatmap.png")
```

### 3.3 `uiuc_propeller_comparison.py`(UIUC 比較)

**目的**: 自設計と UIUC Crazyflie 実測の比較。

**期待する出力**:
1. Crazyflie 45mm vs 自作 36mm のブレード c(R), β(R) 比較プロット
2. Crazyflie 静推力(CT)vs 自作 BEM 予測の比較プロット
3. 表形式の数値比較

**主要な関数構成**:

```python
def load_uiuc_geometry(filepath: str) -> dict:
    """
    Crazyflie geom ファイル(.txt)を読み込み。
    フォーマット:
        r/R   c/R     beta
        0.15  0.2453  29.162
        ...
    Returns: {'r_R': array, 'c_R': array, 'beta_deg': array}
    """
    pass

def load_uiuc_static_data(filepath: str) -> dict:
    """
    UIUC 静推力試験データ読み込み。
    フォーマット:
        RPM        CT        CP
         4026.667  0.139894  0.131208
         ...
    Returns: {'rpm': array, 'CT': array, 'CP': array}
    """
    pass

def ct_to_thrust(ct: float, rpm: float, D: float, rho: float = 1.225) -> float:
    """
    CT から推力を計算。
    T = CT × ρ × n² × D⁴
    n = RPM / 60
    """
    pass

def compare_thrust_curves_with_uiuc(
    bem_result: dict,
    uiuc_data: dict,
    save_path: str
) -> None:
    """
    自作 BEM 予測 vs Crazyflie UIUC 実測の推力比較プロット。
    """
    pass

def compare_geometry_with_uiuc(
    self_geom: BladeGeometry,
    uiuc_geom: dict,
    save_path: str
) -> None:
    """
    自作 vs Crazyflie のブレード形状比較プロット。
    """
    pass

if __name__ == "__main__":
    # 1. Crazyflie データ読み込み
    cf_geom_data = load_uiuc_geometry("uiuc_propeller_data/volume-2/data/cfnq_45p1_geom.txt")
    cf_static_data = load_uiuc_static_data("uiuc_propeller_data/volume-2/data/cfnq_45p1_static_1322rd.txt")
    
    # 2. 自作プロペラのスケーリング後形状
    self_geom = load_blade_geometry_from_crazyflie(D=0.036)
    
    # 3. 形状比較プロット
    compare_geometry_with_uiuc(self_geom, cf_geom_data, "plots/geometry_comparison.png")
    
    # 4. 推力比較
    polars = load_airfoil_polars(self_geom.airfoil_name, re_values=[8000, 12000, 16000, 20000])
    bem_result = bem_analysis(self_geom, 22100, V_CRUISE, polars)
    
    # 5. Crazyflie 実測との照合(D・RPM スケーリング含む)
    print("Crazyflie 45mm @ 19,860 rpm 実測:")
    T_cf = ct_to_thrust(cf_static_data['CT'][-1], cf_static_data['rpm'][-1], 0.045)
    print(f"  T = {T_cf*1000:.1f} mN")
    print(f"自作 36mm @ 22,100 rpm BEM 予測: T = {bem_result['thrust']*1000:.1f} mN")
```

---

## 4. Crazyflie UIUC データの活用方法

### 4.1 利用可能なファイル(Nobu が事前取得済み、`drone-projects/uiuc_propeller_data/volume-2/data/`)

**主参考(Crazyflie 45mm)**:
- `cfnq_45p1_geom.txt` ← **最重要(主参考)**
- `cfnq_45p1_static_1322rd.txt` ← **最重要(静推力)**
- `cfnq_45p2_geom.txt`(変種、参考)
- `cfnq_45p2_static_1323rd.txt`(変種、参考)
- `cfnq_45t1_geom.txt`(変種、参考)
- `cfnq_45t1_static_1320rd.txt`(変種、参考)
- `cfnq_45t2_geom.txt`(変種、参考)
- `cfnq_45t2_static_1321rd.txt`(変種、参考)

**補助参考**:
- `pl_57x20_*`(Plantraco 57mm、P/D=0.35、低ピッチ系)
- `gwsdd_2.5x0.8_*`(GWS 63.5mm、低ピッチ)
- `gwsdd_2.5x1_*`(GWS 63.5mm、P/D=0.4)

### 4.2 ファイル形式

**geom ファイル**(例: `cfnq_45p1_geom.txt`):
```
r/R   c/R     beta
0.15  0.2453  29.162
0.20  0.2631  26.789
...
1.00  0.0465  11.537
```
- 1 行目: ヘッダ
- 2 行目以降: 18 行(r/R = 0.15 から 1.00 まで 0.05 刻み)
- スペース区切り

**static データファイル**(例: `cfnq_45p1_static_1322rd.txt`):
```
RPM        CT        CP
 4026.667  0.139894  0.131208
 5046.667  0.149527  0.116677
...
19860.000  0.164534  0.104108
```
- 1 行目: ヘッダ
- 2 行目以降: RPM ごとの静推力測定値
- 注: 静推力 = advance ratio J=0(機体速度ゼロ)での測定

### 4.3 推力換算式

```python
# CT から推力 T への換算
T = CT × ρ × n² × D⁴
# ρ = 1.225 kg/m³
# n = RPM / 60(rev/s)
# D = 0.045 m(Crazyflie の直径)
```

Crazyflie 45mm @ 19,860 rpm の場合:
- CT = 0.165
- n = 331 rps
- T = 0.165 × 1.225 × 331² × 0.045⁴ = 0.0164 N(1.67 g重)

### 4.4 自作プロペラへのスケーリング論理

**形状(c/R, β)はそのまま流用**:
- c/R は無次元なので、D を変えるだけで適用可能
- 例: 自作 D=36mm の 75%R 弦長 = 0.379 × 18mm = 6.82 mm

**推力スケーリング**:
- T ∝ ρ × n² × D⁴(同じ形状の場合)
- ただし **動作 Re が変わるので CT も変わる可能性**
- → BEM で直接計算するのが正解(単純スケーリングは検算用)

---

## 5. 出力する数値・図表

### 5.1 数値出力(コンソール printまたは結果ファイル `.txt`)

各スクリプトの実行で、以下の数値が出力されるようにする:

**`propeller_design_analysis.py`**:
```
=== プロペラ性能解析 ===
形状: D=36mm, P=27mm, P/D=0.75
動作Re域(@22.1K rpm): 8,200 - 18,500
巡航推力(22.1K rpm, 7.5m/s): __.__ mN (目標 26.5 mN)
巡航効率: 0.___
最大効率時の前進比 J: 0.___
最大推力(40K rpm, 0m/s 静止): ___ mN
```

**`propeller_structural_analysis.py`**:
```
=== 強度解析 ===
ブレード遠心力(40K rpm 最悪): __._ N
根本断面: 5.0 × 1.5 mm² = 7.5 mm²

40K最悪: σ_T=__.__ MPa, σ_B=_.__ MPa, σ_total=__.__ MPa, SF=_._
27K最大動作: σ_total=__.__ MPa, SF=_._
22K巡航: σ_total=__.__ MPa, SF=__._

翼端マッハ数(40K rpm): 0.___
軸-ハブ嵌合 SF: ___
```

**`uiuc_propeller_comparison.py`**:
```
=== UIUC 比較 ===
Crazyflie 45mm @ 19,860 rpm: T = 16.4 mN(静推力)
自作 36mm @ 22,100 rpm BEM 予測: T = __._ mN(V=7.5m/s 巡航)
推力スケーリング検算: T = __._ mN(静推力スケーリング、参考)
形状一致度: 75%R β差 = _._° (Crazyflie 16.1° vs 自作 17°)
```

### 5.2 図表(PNG 形式、`plots/` に保存)

| # | ファイル名 | 内容 |
|---|---|---|
| 1 | `propeller_3view.png` | プロペラ三面図(正面・側面・上面) |
| 2 | `blade_chord_distribution.png` | c(R) 比較(自作 vs Crazyflie) |
| 3 | `blade_twist_distribution.png` | β(R) 比較(自作 vs Crazyflie) |
| 4 | `airfoil_polars.png` | ブレード断面翼型ポーラ(NeuralFoil) |
| 5 | `thrust_vs_rpm_curves.png` | 推力 vs RPM(V=0,5,7.5,9 m/s の4曲線) |
| 6 | `efficiency_vs_J.png` | η vs J 曲線 |
| 7 | `stress_distribution.png` | ブレード根応力分布 |
| 8 | `sf_sensitivity_heatmap.png` | SF 感度ヒートマップ(m × RPM) |
| 9 | `hub_attachment_detail.png` | ハブ取付詳細図 |
| 10 | `crazyflie_comparison_geom.png` | Crazyflie vs 自作 形状比較 |
| 11 | `crazyflie_comparison_thrust.png` | Crazyflie vs 自作 推力比較 |

---

## 6. ブランチ運用ルール(厳守)

🚨 **重要**: Nobu さんは別のリポジトリも管理しているため、**Claude Code が勝手にブランチを切らない** ようにすること。

### 6.1 やってよいこと

- ✅ 既存のブランチ(`main` or 現在のブランチ)で新規ファイルを作成
- ✅ 新規ファイルのコミット(Nobu さんの確認後)

### 6.2 やってはいけないこと

- ❌ **新規ブランチの作成**(`git checkout -b ...`)
- ❌ **既存ファイルの大幅な変更**(主翼・尾翼の解析スクリプトに手を加えない)
- ❌ **Nobu の確認なしの push**

### 6.3 ブランチが必要な場合

「ブランチを切ったほうが良いか?」と Nobu さんに確認する。指示があるまで現在のブランチで作業。

---

## 7. コーディング規約

### 7.1 関数の書き方

- **docstring 必須**: 関数の目的・引数・戻り値を記述
- **型ヒント推奨**: `def func(x: float) -> float:` の形式
- **コメント**: 「なぜそうするか」を書く(「何をしているか」はコードで分かる)

例:
```python
def compute_centrifugal_force(blade_mass: float, r_cg: float, omega: float) -> float:
    """
    ブレード遠心力を計算。
    
    プロペラ設計の支配荷重(推力の数百倍になる)。
    根本断面の引張応力評価に使用。
    
    Args:
        blade_mass: ブレード質量 [kg]
        r_cg: ブレード重心位置 [m](通常 0.75 × R)
        omega: 角速度 [rad/s]
    
    Returns:
        遠心力 F [N]
    
    例:
        >>> compute_centrifugal_force(0.0003, 0.0135, 4189)
        71.0
    """
    return blade_mass * r_cg * omega**2
```

### 7.2 設計値の管理

ハードコーディングは避け、`config` または `dataclass` で管理:

```python
# 良い例
from dataclasses import dataclass

@dataclass
class DesignConstants:
    rpm_max: float = 40000
    diameter: float = 0.036
    # ...

# 悪い例(数値が散らばる)
def some_function():
    rpm = 40000  # ← どこで定義したか不明
    D = 0.036    # ← 後で変更しづらい
```

### 7.3 主翼スクリプトとの統一感

`mainwing_structural_analysis.py` を参考にして、以下のスタイルを継承:
- 関数構成
- 出力形式(printの書式)
- プロット様式(matplotlib のスタイル設定)
- ファイル冒頭のモジュール docstring

### 7.4 numpy 注意

- numpy 2.x の `np.trapezoid` を使用(`np.trapz` は削除済み)

### 7.5 プロット保存先

```python
import os
os.makedirs("plots", exist_ok=True)
plt.savefig("plots/propeller_3view.png", dpi=150, bbox_inches='tight')
```

PNG は `plots/` に保存。後で Nobu が手動で Google Drive にコピー(`ドローン/data関連/グラフ/2026-05-15_propeller_design/`)。

---

## 8. 自己レビュー観点(実装完了時にチェック)

### 8.1 数値の整合性

- [ ] 巡航推力 BEM 予測 = 必要推力 26.5 mN ± 30% 以内(理想)
- [ ] SF(40K rpm) ≥ 3.5、SF(27K rpm) ≥ 6 を維持
- [ ] 翼端マッハ数 < 0.3(40K rpm 時)
- [ ] 軸-ハブ嵌合 SF ≥ 100

### 8.2 スクリプト品質

- [ ] 各関数に docstring が記述されている
- [ ] 設計値が `dataclass` または `config` で管理されている(ハードコーディングなし)
- [ ] エラー処理が入っている(ファイル読み込み失敗、収束失敗など)
- [ ] 実行コマンド一発で全プロット生成されること(`python propeller_design_analysis.py`)

### 8.3 UIUC データ活用

- [ ] Crazyflie geom データを正しく読み込めている
- [ ] Crazyflie 静推力との比較プロットが生成されている
- [ ] スケーリングによる推力予測が計算書に書ける形になっている

### 8.4 図表の品質

- [ ] 全 11 図表が `plots/` に生成されている
- [ ] 軸ラベル・タイトル・凡例が日本語または英語で適切
- [ ] 解像度 150 dpi 以上で出力されている
- [ ] 計算書(Markdown)に貼り付けても文字が潰れない

### 8.5 設計書との整合性

- [ ] 設計書 v0.3 §3 の凍結値(D=36, P=27 など)が変更されていない
- [ ] 設計書 v0.3 §5 の SF 数値(3.9 / 8.9 / 12.9)と一致している
- [ ] 設計書 v0.3 §4.7 の推力スケーリング検算(0.83 g重)と整合

### 8.6 完了したら

完了したら以下を Nobu さんに報告:
1. 生成したファイル一覧(.py × 3、.png × 11)
2. 数値出力の主要結果(巡航推力、各 SF など)
3. 想定外の知見・残課題(あれば)

その後、計算書 §5 原稿(.md)作成フェーズへ。

---

## 9. 残課題・既知の制限

### 9.1 ブレード質量未測定

`M_BLADE_TENTATIVE = 0.3e-3 # kg` は仮値。
明日 5/14 に Nobu が 0.01g 精度の天秤で測定 → 実測値で更新。
**Claude Code は最初は仮値で実装し、感度解析で 0.2-0.5g 範囲を評価しておく**(後で実測値が入ったら数値だけ差し替えれば良い設計にする)。

### 9.2 翼型未確定

NeuralFoil で 4 候補(S1223, E205, 薄キャンバ平板, AG12)を動作Re域でスクリーニング → 1 つに確定。
**Claude Code が翼型選定まで担当**。ただし最終判断は Nobu さんと共有(分かりやすい比較プロットを提供)。

### 9.3 BEM 自作実装の精度

XROTOR/QPROP を使わず Python で自作するため、絶対値精度は ±20% 程度の想定。
**UIUC Crazyflie データとの比較で「予測トレンドが実測と整合」を示せれば、絶対値の不確実性は受容**。

### 9.4 静推力 vs 飛行中推力

UIUC Crazyflie データは **静推力(V=0)** のみ。
飛行中(V=7.5 m/s)の推力は静推力より大きくなる(プロペラ理論の一般則)が、その正確な係数は不明。
**BEM で V>0 の推力を直接計算** することで対応。

---

## 10. 質問・連絡先

- 設計値の変更・追加が必要な場合 → **管理スレッド(本会話)に確認**
- 物理的に解析を進められない不確実性が出た場合 → **管理スレッドに報告**
- Nobu さんの作業(天秤確認、追加実測)が必要な場合 → **本指示書に明記して管理スレッドに連絡**

---

**Claude Code への期待**: 設計書 v0.3 の数値根拠を、再現可能な Python スクリプトとして実装し、計算書 §5 の数値・図表を生成する。
**期限**: 2026-05-15(金)夕方。
**品質基準**: 主翼 §2 の `mainwing_structural_analysis.py` と同等の網羅性・docstring の充実度。

頑張ってください!
