# monolis solid

- monolis structural analysis solver

## 入力ファイル

- input.dat
    - 解析パラメータの設定
- node.dat
    - メッシュの節点定義
- elem.dat
    - メッシュの要素定義
- bc.dat
    - ディリクレ境界条件定義
- load.dat
    - ノイマン境界条件定義

## 入力フォーマット詳細

### input.dat

解析パラメータの設定

```
     0  !> isNLGeom (0:off, 1:on)
     5  !> max NR step
200000  !> Young 率
   0.3  !> Poisson 比
   1.0  !> density
```

| 項目 | 説明 | 型 |
| ---- | ---- | ---- |
| isNLGeom | 幾何学的非線形の有効・無効 | 0, 1 (0:off, 1:on) |
| max NR step | Newton-Raphson 反復の上限 | 整数 |
| Young 率 | ヤング率 | 実数 |
| Poisson 比 | ポアソン比 | 実数 |
| density | 密度 | 実数 |

### node.dat

節点の定義

```
{節点数}
{x 座標} {y 座標} {z 座標}
# 節点数分定義
```

### elem.dat

要素の定義

```
{要素数} {要素を構成する節点数}
{節点 id 1} {節点 id 2} {節点 id 3} {節点 id 4} ...
# 要素数分定義
```

### bc.dat

ディリクレ境界条件の定義

```
{境界条件数} {節点の自由度数 (ここでは 3)}
{節点 id} {自由度方向 1, 2, 3} {変位値}
# 境界条件数分定義
```

自由度方向として以下の意味を持つ 1, 2, 3 の値を定義できる。

- `1`: x 自由度
- `2`: y 自由度
- `3`: z 自由度

### load.dat

ノイマン境界条件の定義

```
{境界条件数} {節点の自由度数 (ここでは 3)}
{節点 id} {自由度方向 1, 2, 3} {節点力値}
# 境界条件数分定義
```

自由度方向として以下の意味を持つ 1, 2, 3 の値を定義できる。

- `1`: x 自由度
- `2`: y 自由度
- `3`: z 自由度

