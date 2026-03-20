# `Control law` 論文對照整理稿

## 0. 記號表

### 0.1 時間、軌跡、量測

```text
T_s                : sampling time
\lambda_c          : control pole
k                  : discrete-time index

\mathbf{p}[k]      = [p_x[k], p_y[k], p_z[k]]^T
\mathbf{p}_d[k]    = [x_d[k], y_d[k], z_d[k]]^T

d_x[k]             = x_d[k] - p_x[k]
d_y[k]             = y_d[k] - p_y[k]
d_z[k]             = z_d[k] - p_z[k]

d_{x,m}[k]         : measured x-axis tracking error
d_{z,m}[k]         : measured z-axis tracking error
\bar d_{x,m}[k]    : deterministic component of d_{x,m}[k]
\bar d_{z,m}[k]    : deterministic component of d_{z,m}[k]
d_{x,r}[k]         : random component of d_{x,m}[k]
d_{z,r}[k]         : random component of d_{z,m}[k]
```

### 0.2 幾何、阻力、motion gain

```text
\theta, \phi, p_w  : wall orientation / offset parameters
R                  : probe radius
h[k]               : probe-to-wall distance
\bar h[k]          = h[k] / R

C_{\parallel}(\bar h), C_{\perp}(\bar h)
                   : near-wall correction factors

\gamma_N           : nominal drag coefficient
\gamma_x[k]        = \gamma_N C_{\parallel}(\bar h[k])
\gamma_y[k]        = \gamma_N C_{\parallel}(\bar h[k])
\gamma_z[k]        = \gamma_N C_{\perp}(\bar h[k])

a_x^{(g)}[k]       = T_s / \gamma_x[k]
a_y^{(g)}[k]       = T_s / \gamma_y[k]
a_z^{(g)}[k]       = T_s / \gamma_z[k]

a_{x,m}[k]         : measured x-axis motion gain from variance inversion
a_{z,m}[k]         : measured z-axis motion gain from variance inversion
```

### 0.3 控制力、估測狀態、協方差

```text
\mathbf{f}_d[k]    = [f_{d,x}[k], f_{d,y}[k], f_{d,z}[k]]^T
\mathbf{f}_T[k]    : thermal force

\hat d_{x,1}[k], \hat d_{x,2}[k], \hat d_{x,3}[k]
                   : x-axis delayed error states
\hat x_D[k], \hat{\dot x}_D[k]
                   : x-axis disturbance states
\hat a_x[k], \hat{\dot a}_x[k]
                   : x-axis motion-gain states

\hat d_{z,1}[k], \hat d_{z,2}[k], \hat d_{z,3}[k]
                   : z-axis delayed error states
\hat z_{D,1}[k], \hat z_{D,2}[k]
                   : z-axis disturbance pair
\hat a_{z,1}[k], \hat a_{z,2}[k]
                   : z-axis motion-gain pair

L_x[k], L_z[k]     : estimator feedback matrices
P_{f,x}[k], P_{f,z}[k]
                   : forecast covariance matrices
P_x[k], P_z[k]     : updated covariance matrices
H                  : output matrix
F_{e,x}[k], F_{e,z}[k]
                   : estimator error-state transition matrices
Q_x[k], Q_z[k], R_x[k], R_z[k]
                   : process / measurement covariance matrices
```

---

## 1. 論文主線 A — Discrete-Time Control Law `(6)~(8)`

### 1.1 Paper `(6)`：控制力骨架

```text
f_{d,x}[k]
= \hat a_x^{-1}[k]
  \left(
    x_d[k+1] - x_d[k] + (1-\lambda_c)\hat d_x[k] - \hat x_D[k]
  \right)
```

```text
d_x[k+1] = \lambda_c d_x[k] - e[k]                                      ...(7)

e[k] = a_x[k] f_T[k]
     - (1-\lambda_c)e_{d_x}[k]
     + e_{a_x}[k] f_{d,x}[k]
     + e_{x_D}[k]                                                       ...(8)
```

### 1.2 Chart 現況主線

```text
f_{d,x}[k]
= \left(a_x^{(g)}[k]\right)^{-1}
  \left(
    x_d[k] - x_d[k-1] + (1-\lambda_c) d_x[k-2]
  \right)
```

```text
f_{d,y}[k]
= \left(a_y^{(g)}[k]\right)^{-1}
  \left(
    y_d[k] - y_d[k-1] + (1-\lambda_c) d_y[k-2]
  \right)
```

```text
f_{d,z}[k]
= \left(a_z^{(g)}[k]\right)^{-1}
  \left(
    z_d[k] - z_d[k-1] + (1-\lambda_c)\hat d_{z,3}[k] - \hat z_{D,1}[k-1]
  \right)
```

差異註記：`x/y` 軸目前使用 gain-based 簡化控制；`z` 軸保留 estimator-based control，且分母使用幾何近壁模型得到的 `a_z^{(g)}[k]`，不是 `\hat a_{z,1}[k]`。

### 1.3 x/y/z 誤差定義

```text
d_x[k-2] = x_d[k-2] - p_x[k]
d_y[k-2] = y_d[k-2] - p_y[k]
d_z[k-2] = z_d[k-2] - p_z[k]
```

差異註記：chart 直接以 two-step delayed desired position 搭配當前 noisy measurement 組成誤差；`z` 軸 estimator 注入時使用的是 `d_z[k-2]`。

---

## 2. 論文主線 B — Process Outputs / IIR / Variance-to-Gain `(9)~(13)`

### 2.1 Paper `(9)~(10)`：IIR 與 variance

```text
\bar d_{x,m}[k+1]
= a_{\mathrm{var}} d_{x,m}[k] + (1-a_{\mathrm{var}})\bar d_{x,m}[k]     ...(9)

d_{x,r}[k] = d_{x,m}[k] - \bar d_{x,m}[k]

\sigma^2_{d_{x,r}}[k]
= \overline{d_{x,r}^2}[k] - \left(\overline{d_{x,r}}[k]\right)^2         ...(10)
```

### 2.2 Chart：x 軸 IIR 鏈

```text
d_{x,m}[k] = d_x[k-2]

\bar d_{x,m}[k]
= A_1 d_{x,m}[k-1] + (1-A_1)\bar d_{x,m}[k-1]

d_{x,r}[k] = d_{x,m}[k] - \bar d_{x,m}[k]

\bar d_{x,r}[k]
= A_2 d_{x,r}[k] + (1-A_2)\bar d_{x,r}[k-1]

\tilde d_{x,r}[k] = d_{x,r}[k] - \bar d_{x,r}[k]

\overline{\tilde d_{x,r}}[k]
= A_{22} d_{x,r}[k] + (1-A_{22})\overline{\tilde d_{x,r}}[k-1]

\overline{d_{x,r}^2}[k]
= A_3 d_{x,r}^2[k] + (1-A_3)\overline{d_{x,r}^2}[k-1]

\sigma^2_{d_{x,r}}[k]
= \overline{d_{x,r}^2}[k] - \left(\overline{\tilde d_{x,r}}[k]\right)^2
```

差異註記：chart 的 mean / mean-square recursion 使用 `A_2, A_{22}, A_3` 擴充版 IIR，而非 paper 只寫出的最簡單形式。

### 2.3 Chart：z 軸 IIR 鏈

```text
d_{z,m}[k] = d_z[k-2]

\bar d_{z,m}[k]
= A_1 d_{z,m}[k] + (1-A_1)\bar d_{z,m}[k-1]

d_{z,r}[k] = d_{z,m}[k-1] - \bar d_{z,m}[k]

d_{z,r}[k] \leftarrow \mathrm{clip}(d_{z,r}[k], -6\times10^{-2}, 6\times10^{-2})

\bar d_{z,r}[k]
= A_2 d_{z,r}[k] + (1-A_2)\bar d_{z,r}[k-1]

\tilde d_{z,r}[k] = d_{z,r}[k] - \bar d_{z,r}[k]

\overline{\tilde d_{z,r}}[k]
= A_{22} d_{z,r}[k] + (1-A_{22})\overline{\tilde d_{z,r}}[k-1]

\overline{d_{z,r}^2}[k]
= A_3 d_{z,r}^2[k] + (1-A_3)\overline{d_{z,r}^2}[k-1]

\sigma^2_{d_{z,r}}[k]
= \overline{d_{z,r}^2}[k] - \left(\overline{\tilde d_{z,r}}[k]\right)^2

\sigma^2_{d_{z,r}}[k] \leftarrow \max(\sigma^2_{d_{z,r}}[k], 0)
```

差異註記：`z` 軸使用 bugfix 後的 delayed residual `d_{z,m}[k-1] - \bar d_{z,m}[k]`，並加入 clip 與 nonnegative guard。

### 2.4 Paper `(12)~(13)`：variance-to-gain

```text
\sigma^2_{d_{x,r}}[k]
= \left(2 + \frac{1}{1-\lambda_c^2}\right) 4 k_B T a_x[k]
 + \frac{2}{1+\lambda_c}\sigma^2_{n_x}                                  ...(12)

a_{x,m}[k]
= \frac{
      \sigma^2_{d_{x,r}}[k] - \frac{2}{1+\lambda_c}\sigma^2_{n_x}
   }{
      4k_B T \left(2 + \frac{1}{1-\lambda_c^2}\right)
   }                                                                      ...(13)
```

### 2.5 Chart：x/z motion gain recovery

```text
D(\lambda_c) = 4k_B T \left(2 + \frac{1}{1-\lambda_c^2}\right)

N_x[k] = \sigma^2_{d_{x,r}}[k] - \frac{2}{1+\lambda_c}\sigma^2_{n_x}
N_z[k] = \sigma^2_{d_{z,r}}[k] - \frac{2}{1+\lambda_c}\sigma^2_{n_z}

N_x[k] \leftarrow \max(N_x[k], 0)
N_z[k] \leftarrow \max(N_z[k], 0)

a_{x,m}^{\mathrm{dyn}}[k] = \frac{N_x[k]}{D(\lambda_c)} \cdot \mathrm{unit\_scale}
a_{z,m}^{\mathrm{dyn}}[k] = \frac{N_z[k]}{D(\lambda_c)} \cdot \mathrm{unit\_scale}

a_{x,m}[k] = 2.5 \, a_{x,m}^{\mathrm{dyn}}[k]
a_{z,m}[k] = A_m \, a_{z,m}^{\mathrm{dyn}}[k]
```

差異註記：chart 對 `(13)` 加入單位換算、nonnegative guard，以及 `2.5` / `A_m` scaling；`z` 軸 measured gain 經過額外縮放，不等於 paper 的直接反解式。

---

## 3. 論文主線 C — 7-State Estimator / Covariance Recursion `(14)~(21)`

### 3.1 Paper `(14)~(18)`：7-state process 與 estimator

```text
\mathbf{x}[k]
=
\begin{bmatrix}
 d_{x,1}[k] & d_{x,2}[k] & d_{x,3}[k] &
 x_D[k] & \dot x_D[k] &
 a_x[k] & \dot a_x[k]
\end{bmatrix}^T
```

```text
d_{x,1}[k+1] = d_{x,2}[k]
d_{x,2}[k+1] = d_{x,3}[k]
d_{x,3}[k+1] = \lambda_c d_{x,3}[k] + (1-\lambda_c)e_{d_x}[k]
               - e_{x_D}[k] - f_{d,x}[k]e_{a_x}[k] - a_x[k]f_T[k]
x_D[k+1]     = x_D[k] + \dot x_D[k]
\dot x_D[k+1]= \dot x_D[k]
a_x[k+1]     = a_x[k] + \dot a_x[k]
\dot a_x[k+1]= \dot a_x[k]                                                  ...(14)
```

```text
H =
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}                                                               ...(15)
```

```text
\hat{\mathbf{x}}[k+1]
= F_e[k]\hat{\mathbf{x}}[k] + L[k]
\begin{bmatrix}
e_{x_1}[k] \\
e_{a_x}[k]
\end{bmatrix}                                                               ...(16)
```

```text
\mathbf{e}[k+1] = F_e[k]\mathbf{e}[k] - L[k]H\mathbf{e}[k] + \mathbf{q}[k] ...(17)
```

### 3.2 Chart：x 軸 estimator

```text
\mathbf{e}_x[k]
=
\begin{bmatrix}
\bar d_{x,m}[k] - \hat d_{x,1}[k-1] \\
a_{x,m}[k] - \hat a_x[k-1]
\end{bmatrix}
```

```text
H =
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}
```

```text
F_{e,x}[k] =
\begin{bmatrix}
0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & -1 & 0 & -f_{d,x}[k-1] & 0 \\
0 & 0 & 0 & 1 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}
```

```text
L_x[k] = P_{f,x}[k-1]H^T \left(HP_{f,x}[k-1]H^T + R_x[k]\right)^{-1}

P_x[k] = (I - L_x[k]H)P_{f,x}[k-1]

P_{f,x}[k] = F_{e,x}[k] P_x[k] F_{e,x}^T[k] + Q_x[k]
```

```text
\hat d_{x,1}[k] = \hat d_{x,2}[k-1] + \mathrm{inj}_1[k]
\hat d_{x,2}[k] = \hat d_{x,3}[k-1] + \mathrm{inj}_2[k]
\hat d_{x,3}[k] = \lambda_c \hat d_{x,3}[k-1] + \mathrm{inj}_3[k]

\hat x_D[k]      = \hat x_D[k-1] + \hat{\dot x}_D[k-1] + \mathrm{inj}_4[k]
\hat{\dot x}_D[k]= \hat{\dot x}_D[k-1] + \mathrm{inj}_5[k]

\hat a_x[k]      = \hat a_x[k-1] + \hat{\dot a}_x[k-1] + \mathrm{inj}_6[k]
\hat{\dot a}_x[k]= \hat{\dot a}_x[k-1] + \mathrm{inj}_7[k]
```

差異註記：x 軸 estimator 基本上延續 paper `(14)~(21)`；目前 `x` 控制器輸出未直接使用 `\hat a_x[k]` / `\hat x_D[k]`，而是改用 near-wall 幾何 gain。

### 3.3 Chart：z 軸 estimator

```text
\mathbf{e}_z[k]
=
\begin{bmatrix}
d_z[k-2] - \hat d_{z,1}[k-1] \\
a_{z,m}[k-1] - \hat a_{z,1}[k-1]
\end{bmatrix}
```

```text
H =
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}
```

```text
F_{e,z}[k] =
\begin{bmatrix}
0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & -1 & 0 & -f_{d,z}[k-1] & 0 \\
0 & 0 & 0 & 1+\beta & -\beta & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1+\beta & -\beta \\
0 & 0 & 0 & 0 & 0 & 1 & 0
\end{bmatrix}
```

```text
R_z[k] = \mathrm{diag}(r_{11}[k], r_{22}[k])
Q_z[k] = \mathrm{diag}(q_{11}[k], q_{22}[k], q_{33}[k], q_{44}[k], q_{55}[k], q_{66}[k], q_{77}[k])
```

```text
L_z[k] = P_{f,z}[k-1]H^T \left(HP_{f,z}[k-1]H^T + R_z[k]\right)^{-1}     ...(19)

P_z[k] = \lambda_F^{-1}(I - L_z[k]H)P_{f,z}[k-1]                          ...(20')

P_{f,z}[k] = F_{e,z}[k] P_z[k] F_{e,z}^T[k] + Q_z[k]                      ...(21)
```

```text
\hat d_{z,1}[k] = \hat d_{z,2}[k-1] + \mathrm{inj}_1[k]
\hat d_{z,2}[k] = \hat d_{z,3}[k-1] + \mathrm{inj}_2[k]
\hat d_{z,3}[k] = \lambda_c \hat d_{z,3}[k-1] + \mathrm{inj}_3[k]

\hat z_{D,1}[k] = (1+\beta)\hat z_{D,1}[k-1] - \beta \hat z_{D,2}[k-1] + \mathrm{inj}_4[k]
\hat z_{D,2}[k] = \hat z_{D,1}[k-1] + \mathrm{inj}_5[k]

\hat a_{z,1}[k] = (1+\beta)\hat a_{z,1}[k-1] - \beta \hat a_{z,2}[k-1] + \mathrm{inj}_6[k]
\hat a_{z,2}[k] = \hat a_{z,1}[k-1] + \mathrm{inj}_7[k]
```

差異註記：`z` 軸將 paper 的兩個 second-order random-walk pair 改寫成 `\beta` 耦合形式；此外加入 forgetting factor `\lambda_F`。

### 3.4 Chart：數值保護

```text
S[k] = H P_{f,\bullet}[k-1] H^T + R_{\bullet}[k]
S[k] \leftarrow \tfrac12(S[k] + S^T[k]) + \varepsilon I

P_{\bullet}[k]   \leftarrow \tfrac12(P_{\bullet}[k] + P_{\bullet}^T[k])
P_{f,\bullet}[k] \leftarrow \tfrac12(P_{f,\bullet}[k] + P_{f,\bullet}^T[k])
```

差異註記：paper `(19)~(21)` 沒有寫出正則化與 symmetry enforcement；chart 實作中明確加入。

---

## 4. Chart-only 前後處理

### 4.1 初始化

```text
\hat d_{x,1}[0] = \hat d_{x,2}[0] = \hat d_{x,3}[0] = 0
\hat d_{z,1}[0] = \hat d_{z,2}[0] = \hat d_{z,3}[0] = 0
\hat x_D[0] = \hat{\dot x}_D[0] = 0
\hat z_{D,1}[0] = \hat z_{D,2}[0] = 0
\hat a_x[0] = 0.012
\hat a_{z,1}[0] = 0.012/(1+\beta)
P_{f,x}[0], P_{f,z}[0] = \text{diagonal initial covariances}
```

差異註記：這些初值屬於 chart implementation，paper 主文未固定指定。

### 4.2 Near-wall 幾何與 thermal-force 中介量

```text
\mathbf{\hat w}
=
\begin{bmatrix}
\cos\theta\sin\phi \\
\sin\theta\sin\phi \\
\cos\phi
\end{bmatrix}

h[k] = \mathbf{p}^T[k]\mathbf{\hat w} - p_w
\bar h[k] = h[k] / R

C_{\parallel}(\bar h[k]), C_{\perp}(\bar h[k])
\rightarrow
\gamma_x[k], \gamma_y[k], \gamma_z[k]
\rightarrow
a_x^{(g)}[k], a_y^{(g)}[k], a_z^{(g)}[k]
```

```text
\mathbf{f}_T[k] \sim \mathcal{N}\!\left(0, \Sigma_T[k]\right)
```

差異註記：此段屬於 3D near-wall plant-side augmentation，非 paper `(6)~(21)` 估測主線的一部分，但決定了 chart 中 gain 與 thermal statistics 的來源。

### 4.3 量測雜訊注入

```text
p_x[k] \leftarrow p_x[k] + n_x[k], \quad n_x[k]\sim\mathcal{N}(0,\sigma_{n_x}^2)
p_y[k] \leftarrow p_y[k] + n_y[k], \quad n_y[k]\sim\mathcal{N}(0,\sigma_{n_x}^2)
p_z[k] \leftarrow p_z[k] + n_z[k], \quad n_z[k]\sim\mathcal{N}(0,\sigma_{n_z}^2)
```

差異註記：paper 只使用 measurement-noise variance 進入 `(12)` 與 `(19)`，chart 則直接在量測端顯式注入 noise。

### 4.4 State update

```text
x_d[k-2] \leftarrow x_d[k-1], \quad x_d[k-1] \leftarrow x_d[k]
y_d[k-2] \leftarrow y_d[k-1], \quad y_d[k-1] \leftarrow y_d[k]
z_d[k-2] \leftarrow z_d[k-1], \quad z_d[k-1] \leftarrow z_d[k]
```

```text
\hat{\mathbf{x}}_x[k-1] \leftarrow \hat{\mathbf{x}}_x[k], \quad
\hat{\mathbf{x}}_z[k-1] \leftarrow \hat{\mathbf{x}}_z[k]
```

```text
P_{f,x}[k-1] \leftarrow P_{f,x}[k], \quad
P_{f,z}[k-1] \leftarrow P_{f,z}[k]
```

```text
f_{d,x}[k-1] \leftarrow f_{d,x}[k], \quad
f_{d,z}[k-1] \leftarrow f_{d,z}[k]

d_{x,m}[k-1] \leftarrow d_{x,m}[k], \quad
d_{z,m}[k-1] \leftarrow d_{z,m}[k], \quad
a_{z,m}[k-1] \leftarrow a_{z,m}[k]
```

差異註記：chart 採 explicit persistent-state update，對應實際 FPGA / chart 逐步執行順序。

---

## 5. 符號對照表

| 數學符號 | chart 原始變數 |
| --- | --- |
| `\lambda_c` | `lamdaC` |
| `T_s` | `Ts` |
| `p_x, p_y, p_z` | `px_k`, `py_k`, `pz_k` |
| `x_d, y_d, z_d` | `xd_k`, `yd_k`, `zd_k` |
| `d_x[k-2], d_y[k-2], d_z[k-2]` | `dx_k2`, `dy_k2`, `dz_k2` |
| `d_{x,m}` | `dxm_k` |
| `d_{z,m}` | `dzm_k` |
| `\bar d_{x,m}` | `dxm_bar_k` |
| `\bar d_{z,m}` | `dzm_bar_k` |
| `d_{x,r}` | `dxr_k` |
| `d_{z,r}` | `dzr_k` |
| `\bar d_{x,r}` | `dxr_bar_k` |
| `\bar d_{z,r}` | `dzr_bar_k` |
| `\overline{\tilde d_{x,r}}` | `dxrr_bar_k` |
| `\overline{\tilde d_{z,r}}` | `dzrr_bar_k` |
| `\overline{d_{x,r}^2}` | `dxr2_bar_k` |
| `\overline{d_{z,r}^2}` | `dzr2_bar_k` |
| `\sigma^2_{d_{x,r}}` | `sigma2_dxr` |
| `\sigma^2_{d_{z,r}}` | `sigma2_dzr` |
| `a_{x,m}` | `axm_k` |
| `a_{z,m}` | `azm_k` |
| `a_x^{(g)}` | `mgain_x` |
| `a_y^{(g)}` | `mgain_y` |
| `a_z^{(g)}` | `mgain_z` |
| `f_{d,x}, f_{d,y}, f_{d,z}` | `fd_X`, `fd_Y`, `fd_Z` |
| `\hat d_{x,1}, \hat d_{x,2}, \hat d_{x,3}` | `dx1_hat_k`, `dx2_hat_k`, `dx3_hat_k` |
| `\hat d_{z,1}, \hat d_{z,2}, \hat d_{z,3}` | `dz1_hat_k`, `dz2_hat_k`, `dz3_hat_k` |
| `\hat x_D, \hat{\dot x}_D` | `xD_hat_k`, `dxD_hat_k` |
| `\hat z_{D,1}, \hat z_{D,2}` | `zD1_hat_k`, `zD2_hat_k` |
| `\hat a_x, \hat{\dot a}_x` | `ax_hat_k`, `dax_hat_k` |
| `\hat a_{z,1}, \hat a_{z,2}` | `az1_hat_k`, `az2_hat_k` |
| `L_x, L_z` | `Lx`, `Lz` |
| `P_{f,x}, P_{f,z}` | `Pfx_k`, `Pfz_k` |
| `P_x, P_z` | `Px`, `Pz` |
| `F_{e,x}, F_{e,z}` | `Fe` in x block / `Fe` in z block |
| `Q_x, Q_z` | `Q` in x block / `Q` in z block |
| `R_x, R_z` | `RR` in x block / `RR` in z block |
| `\theta, \phi, p_w` | `theta`, `phi`, `pz` |
| `R` | `R` |
| `C_{\parallel}, C_{\perp}` | `c_par`, `c_perp` |
| `\gamma_x, \gamma_y, \gamma_z` | `gamma_x`, `gamma_y`, `gamma_z` |

---

## 6. paper-vs-chart 差異表

| paper 對應式號 | chart 實作 | 差異類型 | 是否影響主線 |
| --- | --- | --- | --- |
| `(6)` | `x/y` 軸使用 `a_x^{(g)}, a_y^{(g)}` 與 delayed error 的簡化控制律 | x/y 簡化 | 是 |
| `(6)` | `z` 軸分母使用 `mgain_z`，不是 `az1_hat_k` | 控制律替換 | 是 |
| `(7)(8)` | `z` 軸誤差注入以 `dz_k2` 為量測主入口 | delayed measurement 實作化 | 是 |
| `(9)(10)` | IIR 拆成 `Avar, Avar2, Avar22, Avar3` 四組係數 | 參數化延伸 | 是 |
| `(9)(10)` | `z` 軸 `dzm_bar_k` / `dzr_k` 使用 bugfix 後公式 | bugfix | 是 |
| `(9)(10)` | `z` 軸 residual 加入 clip，variance 加入 `max(\cdot,0)` | 數值保護 | 是 |
| `(12)(13)` | 反推 gain 後再乘 `2.5` 或 `Am_scaling` | scaling 延伸 | 是 |
| `(13)` | `x/z` measured gain 需要顯式單位換算 | 單位實作 | 是 |
| `(14)(18)` | `z` 軸的 disturbance / gain 狀態改為 `\beta` 耦合 second-order pair | 狀態模型擴充 | 是 |
| `(15)` | `H` 與 paper 相同，但 `z` 軸第二量測使用 `azm_k1` | delayed output 實作化 | 是 |
| `(19)~(21)` | `z` 軸加入 forgetting factor `lamdaF` | 濾波器延伸 | 是 |
| `(19)~(21)` | `S`, `P`, `Pf` 都做 symmetry enforcement 與 `epsI` regularization | 數值保護 | 是 |
| `none` | chart 先計算 near-wall 幾何、`c_par/c_perp`、`mgain_x/y/z` | 3D 擴充 | 是 |
| `none` | chart 顯式產生 thermal force `ft_3D` 與 measurement noise | plant/measurement 注入 | 否 |
| `none` | x 軸 estimator 還在跑，但最終控制未完整採用 `\hat a_x, \hat x_D` | 部分保留支線 | 否 |
