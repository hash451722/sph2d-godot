# GodotEngineによるSPH法を使った非圧縮性流体表現


## 1. 支配方程式 (Governing equation)

**連続の式 (Equation of continuity)**

$$
\frac{D \rho}{D t} + \rho \nabla \cdot {\bf v} = 0 \tag{1-1}
$$

**ナビエ・ストークス方程式 (Navier-Stokes equation)**

$$
\frac{D {\bf v}}{D t} = -\frac{1}{\rho} \nabla p + \frac{\mu}{\rho} \nabla^2 {\bf v} + {\bf f} \tag{1-2}
$$

ここで，$\rho$は密度，${\bf v}$は速度ベクトル，$p$は圧力，$\mu$は粘性係数，${\bf f}$は外力による加速度ベクトルを表す．左辺の$\frac{D}{D t}$は移動する流体粒子がもつ物理量の時間変化率のことでラグランジュ微分と呼ばれる．

$$
\frac{D}{Dt} = \frac{\partial}{\partial t} + {\bf v} \cdot \nabla \tag{1-3}
$$


## 2. 物理量の離散化

### 2.1 物理量の積分表現
SPH法では物理量$\phi$は，位置$\bf r$の物理量$\phi$をその周辺の物理量に重み関数$W$を掛けて体積積分したものとして表す．

$$
\phi({\bf r}) = \int \phi({\bf r}') W({\bf r}-{\bf r}') dV \tag{2-1}
$$

ここで，重み関数$W$は距離0 $(r=0)$で最大値となる対象形状をしている．これは物理量$\phi$を計算するとき，距離の近い粒子ほどその粒子の影響を大きく受けることを意味している．さらに、距離$h$より離れると関数はゼロとなり，距離$h$より遠い粒子の影響は全く受けないことを表している．この距離$h$を影響半径と呼ぶ．また，重み関数$W$は積分すると$1$になり、その他にも数学的な性質を有する．このような重み関す$W$を特にカーネル関数(Kernel functions)と呼ぶ．

### 2.2 物理量の離散表現
SPH法では式(2-1)を次式のように離散化する．


$$
\phi({\bf r}_i) = \sum_j \frac{m_j}{\rho_j} \phi({\bf r}_j) W({\bf r}_i-{\bf r}_j) \tag{2-2}
$$

ここで，添字$i$は着目粒子，添字$j$は着目粒子以外の粒子を表す．$m_j$，$\rho_j$は粒子jの質量、密度で表しており、体積$dV_j = m_j / \rho_j$の関係を考慮している．つまり着目粒子$i$を中心とした影響半径$h$内の粒子$j$が粒子$i$の物理量$\phi$を決定する．

### 2.3 カーネル関数
SPH法で用いられるカーネル関数は多数提案されているが，ここではPoly6カーネル関数を使用する．Poly6カーネル関数とその勾配は次のように表される．定数$\alpha$は2次元と3次元の計算で異なる値を用いる．
$$
W_{poly6}({\bf r}) =
    \begin{cases}
        \; \alpha ( h^2 - r^2)^3 & r<h \\
        \; 0 & r \ge h
    \end{cases}
\tag{2-3}
$$
$$
\nabla W_{poly6}({\bf r}) =
    \begin{cases}
        \; -6 \alpha ( h^2 - r^2)^2 {\bf r} & r<h \\
        \; 0 & r \ge h
    \end{cases}
\tag{2-4}
$$
$$
\alpha =
    \begin{cases}
        \; 4 / (\pi h^8) & for \space 2D\\
        \; 315/(64 \pi h^9) & for \space 3D
    \end{cases}
\tag{2-5}
$$


### 2.4 コード
カーネル関数のGDScriptでの実装例を示す．クラスは影響半径$h$と定数$\alpha$をメンバ変数として持つ．カーネル関数を計算するメソッドkernelを定義する．引数は中心からの距離$r$で、(2-3)式で与えられる値を返す．そして、勾配を計算するgradientメソッドを定義する．勾配はベクトルなので、引数には距離ベクトルrvを入力し，(2-4)式で計算されたベクトルを返す．


```GDScript
# Poly6 kernel
class Kernel:
	var h: float  # Radius of influence
	var alpha: float
  
	func _init(h_arg: float) -> void:
		h = h_arg
		alpha = 4.0 / (PI * pow(h, 8))

	func kernel(r: float) -> float:
		if (r < h):
			return alpha * pow(h * h - r * r, 3)
		else:
			return 0.0
  
	func gradient(rv: Vector2) -> Vector2:
		var r: float = rv.length()
		if r < h:
			var c: float = -6.0 * alpha * pow(h * h - r * r, 2)
			return Vector2(c * rv.x, c * rv.y)
		else:
			return Vector2.ZERO
```


## 3. 粒子クラス
粒子を表現するParticleクラスの実装例を示す．


```GDScript
class Particle:
	var position := Vector2.ZERO
	var velocity := Vector2.ZERO
	var acceleration := Vector2.ZERO
	var force := Vector2.ZERO
	var pressure: float = 0.0
	var density: float = 0.0
	var mass: float = 0.0
	var active: bool = true

	func _init(pos := Vector2.ZERO):  # constructor
		position = pos
```


## 4. 密度の計算式
SPH法における密度は式(2-2)で示した物理量$\phi$の離散式に密度$\rho$に代入して計算する．
$$
\begin{split} \rho({\bf r}_i) &= \sum_j \frac{m_j}{\rho_j} \rho({\bf r}_j) W({\bf r}_i-{\bf r}_j) \\ &= \sum_j m_j W({\bf r}_i-{\bf r}_j) \end{split} \tag{4-1}
$$


## 5. 圧力の計算式
圧力$p$は状態方程式から求める．状態方程式は様々提案されているが，ここでは基準密度$\rho_0$との密度差に比例した圧力が発生すると仮定した次式を用いる．
$$
p = k (\rho - \rho_0) \tag{5-1}
$$
ここで，$k$は気体定数である．


## 6. ナビエ・ストークス方程式の離散化
粒子の運動方程式はナビエ・ストークス方程式で表せる．右辺の第一項は圧力項，第二項は粘性項，第三項は外力項である．
$$
\frac{D {\bf v}}{D t} = \underbrace{-\frac{1}{\rho} \nabla p}_{圧力項} + \underbrace{\frac{\mu}{\rho} \nabla^2 {\bf v}}_{粘性項} + \underbrace{{\bf f}}_{外力項} \tag{6-1}
$$

### 6.1 微分項の離散化
SPH法における微分項の離散化は式(2-1)より，物理量$\phi$の勾配は
$$
\nabla \phi({\bf r}) = \int \nabla_{r'} \phi({\bf r}') W({\bf r}-{\bf r}') dV \tag{6-2}
$$

上式を式変形することで，物理量$\phi$の勾配は次式で表される．
$$
\nabla \phi({\bf r}) = \int \phi({\bf r}') \nabla_{r} W({\bf r}-{\bf r}') dV \tag{6-3}
$$

離散形では次式となる．
$$
\nabla \phi({\bf r}_i) = \sum_j \frac{m_j}{\rho_j} \phi({\bf r}_j) \nabla W({\bf r}_i-{\bf r}_j) \tag{6-4}
$$
物理量$\phi$の勾配はカーネル関数$W$の勾配で表現することができる．


### 6.2 圧力項の離散化
式(6-4)を用いて圧力項を離散化すると次式となる．
$$
-\frac{1}{\rho} \nabla p = -\frac{1}{\rho_i} \sum_j \frac{m_j}{\rho_j} p_j \nabla W({\bf r}_i-{\bf r}_j) \tag{6-5}
$$
2個の粒子$i$と$j$を考えた場合，それぞれに圧力$p_i$と$p_j$が異なることになり，双方からみた力が一致しない．これは作用反作用の法則に反することになるので，離散式として2個の粒子に働く力は大きさが同じで向きが逆の反対称性を持つ必要がある．
そこで，次式の関係を使って圧力項を書き換える．

$$
\nabla \left( \frac{p}{\rho} \right) = \frac{\nabla p}{\rho} - \frac{p}{\rho^2} \nabla \rho \tag{6-6}
$$
$$
-\frac{1}{\rho} \nabla p = -\left( \nabla \left( \frac{p}{\rho} \right) + \frac{p}{\rho^2} \nabla \rho \right) \tag{6-7}
$$
右辺を離散化する．
$$
\begin{split} -\frac{1}{\rho} \nabla p &= -\left( \sum_j \frac{m_j}{\rho_j} \left( \frac{p_j}{\rho_j} \right) \nabla W({\bf r}_i-{\bf r}_j) + \frac{p_i}{\rho_i^2} \sum_j \frac{m_j}{\rho_j} \rho_j \nabla W({\bf r}_i-{\bf r}_j) \right)\\ &= - \sum_j m_j \left( \frac{p_j}{\rho_j^2} +  \frac{p_i}{\rho_i^2} \right) \nabla W({\bf r}_i-{\bf r}_j) \end{split} \tag{6-8}
$$

得られた式(6-8)は粒子$i$と$j$の相互に見たとき反対称であり，作用反作用の法則を満たす．


### 6.3　粘性項の離散化
粘性項も圧力項と同様に離散化する．
$$
\frac{\mu}{\rho} \nabla^2 {\bf v} = \sum_j m_j \frac{\mu_i + \mu_j}{\rho_i \rho_j} \frac{({\bf r}_i-{\bf r}_j) \cdot \nabla W({\bf r}_i-{\bf r}_j)}{|{\bf r}_i-{\bf r}_j|^2 + \eta^2} ({\bf v}_i - {\bf v}_j) \tag{6-9}
$$
ここで，$\eta^2 = 0.001h^2$を用いる．
上式も反対称性であり，2階微分がカーネル関数の勾配を用いて表される．


### 6.4 外力項の離散化
外力項は重力加速度を適用する．



## 7. 運動方程式の解法
ナビエ・ストークス方程式の左辺は速度の時間微分の解法として，ここではオイラー陽解法を用いる．式


$$
\frac{\bf u^{k+1}_{i} - u^k_i}{\Delta t} = {\lbrack -\frac{1}{\rho} \nabla p + \nu \nabla^2 \bf u + \bf g \rbrack}^k_i \tag{7.1}
$$
$$
\bf u^{k+1}_i = \bf u^k_i + \Delta t {\lbrack -\frac{1}{\rho} \nabla p + \nu \nabla^2 \bf u + \bf g \rbrack}^k_i \tag{7.2}
$$

$$
\bf r^{k+1}_i = \bf r^k_i + \Delta t \bf u ^{k+1}_i \tag{7.3}
$$

## 近傍粒子の探索



## 参考

[1] Müller, M.; Charypar, D.; Gross, M. Particle-Based Fluid Simulation for Interactive Applications. In Proc. of the 2003 ACM SIGGRAPH/Eurographics Symposium on Computer
Animation. 2003, p. 154–159.


[2] M. Desbrun and M. P. Cani. Smoothed particles: A new
paradigm for animating highly deformable bodies. In
Computer Animation and Simulation ’96 (Proceedings of
EG Workshop on Animation and Simulation), pages 61–76.
Springer-Verlag, Aug 1996.
