# Application of the Finite Difference Method to Solve the Seismic Wave Equation

## Introduction

The finite difference method is a numerical technique used to approximate the solutions of partial differential equations by discretizing them. 
$$
\begin{cases} 
u_{tt} = c^2 (u_{xx} + u_{yy}), & (x, y) \in \Omega, \, t \in (0, T], \\
u(x, y, 0) = I(x, y), & (x, y) \in \Omega, \\
u_t(x, y, 0) = V(x, y), & (x, y) \in \Omega, \\
\frac{\partial u}{\partial n} = 0, & (x, y) \in \partial \Omega, \, t \in (0, T].
\end{cases}
$$

where $ u $ is the displacement, $ c $ is the wave speed, $ \Omega $ is the domain, and $ \partial \Omega $ represents the boundary of the domain.

## Discretization

To solve this PDE using the finite difference method, we discretize the spatial and temporal domains. We create a rectangular mesh over the domain $\Omega = [A, B] \times [C, D]$ with spatial step sizes $\Delta x$ and $\Delta y$, and a temporal step size $\Delta t$.

### Rectangular Mesh

Let $ N_x $ and $ N_y $ be the number of grid points in the $ x $ and $ y $ directions, respectively. The spatial grid points are defined as:

$$ x_i = A + i \Delta x, \quad \Delta x = \frac{|A-B|}{N_x}, \quad i = 0, 1, \ldots, Nx, $$
$$ y_j = C + j \Delta y, \quad \Delta y = \frac{|C-D|}{N_y}, \quad j = 0,1, \ldots, N_y. $$

Similarly, let $ N_t $ be the number of time steps. The temporal grid points are defined as:

$$ t_n = n \Delta t, \quad \Delta t = \frac{T}{N_t}, \quad n = 0,1, \ldots, N_t $$

We denote numerical approximation of $ u(x_i,y_j,t_n) $ as $ u^n_{i,j} $.  

### Finite Difference Approximation

#### Forward Difference

On a computer, derivatives are approximated by finite difference expressions; rearranging gives the forward difference approximation:

$$ \frac{f(x + \Delta x) - f(x)}{\Delta x} = f'(x) + O(\Delta x), $$

where $ O(\Delta x) $ means 'terms of order $\Delta x$', i.e., terms which have size similar to or smaller than $\Delta x$ when $\Delta x$ is small. So the expression on the left approximates the derivative of $ f $ at $ x $, and has an error of size $\Delta x$; the approximation is said to be 'first order accurate'.

#### Backward Difference

Rearranging similarly gives the backward difference approximation:

$$ \frac{f(x) - f(x - \Delta x)}{\Delta x} = f'(x) + O(\Delta x), $$

which is also first order accurate, since the error is of order $\Delta x$.

#### Centered Difference

Combining the forward and backward difference approximations gives the centered difference approximation:

$$ \frac{f(x + \Delta x) - f(x - \Delta x)}{2\Delta x} = f'(x) + O(\Delta x^2), $$

which is 'second order accurate', because the error this time is of order $\Delta x^2$.

#### Second Derivative Centered Difference

Adding the forward and backward difference expressions and using Taylor's formula gives:

$$ f(x + \Delta x) + f(x - \Delta x) = 2f(x) + \Delta x^2 f''(x) + \frac{\Delta x^4}{12} f^{(4)}(x) + \cdots $$

Rearranging this therefore gives the centered difference approximation to the second derivative:

$$ \frac{f(x + \Delta x) - 2f(x) + f(x - \Delta x)}{\Delta x^2} = f''(x) + O(\Delta x^2), $$

which is second order accurate.

Using these finite difference approximations, we can discretize the seismic wave equation.

### Discretized Wave Equation

The second-order central difference approximations for the second derivatives are:

$$ u_{xx} \approx \frac{u^n_{i+1,j} - 2u^n_{i,j} + u^n_{i-1,j}}{(\Delta x)^2}, $$ 
$$ u_{yy} \approx \frac{u^n_{i,j+1} - 2u^n_{i,j} + u^n_{i,j-1}}{(\Delta y)^2}, $$

$$ u_{tt} \approx \frac{u_{i,j}^{n+1} - 2u_{i,j}^n + u_{i,j}^{n-1}}{(\Delta t)^2}. $$

Substituting these approximations into the wave equation, we get the finite difference scheme:

$$\frac{u_{i,j}^{n+1} - 2u_{i,j}^n + u_{i,j}^{n-1}}{(\Delta t)^2} = c^2\frac{u_{i+1,j}^{n} - 2u_{i,j}^n + u_{i-1,j}^{n}}{(\Delta x)^2} + c^2\frac{u_{i,j+1}^{n} - 2u_{i,j}^n + u_{i,j-1}^{n}}{(\Delta y)^2} $$

$$ u_{i,j}^{n+1} = 2u_{i,j}^n - u_{i,j}^{n-1} + \left( \frac{c \Delta t}{\Delta x} \right)^2 (u_{i+1,j}^n - 2u_{i,j}^n + u_{i-1,j}^n) + \left( \frac{c \Delta t}{\Delta y} \right)^2 (u_{i,j+1}^n - 2u_{i,j}^n + u_{i,j-1}^n). \tag{**} $$


### Initial Conditions

The initial conditions are specified as:

$$ I(x_i, y_j)=u(x_i, y_j, 0) = u^0_{i,j},  $$
$$ V(x_i, y_j) = u_t(x_i, y_j, 0) = \frac{u_{i, j}^1 - u_{i, j}^{-1}}{2\Delta t} \; \Longrightarrow \; u_{i,j}^{-1} = u_{i,j}^1 - 2\Delta t V(x_i, y_j). $$

To get rid of 'fictional' value $ u_{i,j}^{-1} $ we use equation $ (**) $ with $ n=0 $ and obtain

$$ u_{i,j}^{1} = 2u_{i,j}^0 - u_{i,j}^1 + 2\Delta t V(x_i, y_j) + \left( \frac{c \Delta t}{\Delta x} \right)^2 (u_{i+1,j}^0 - 2u_{i,j}^0 + u_{i-1,j}^0) + \left( \frac{c \Delta t}{\Delta y} \right)^2 (u_{i,j+1}^0 - 2u_{i,j}^0 + u_{i,j-1}^0), $$

$$u_{i,j}^{1} = u_{i, j}^0  + \Delta t V(x_i, y_j) + \frac{1}{2}\left( \frac{c \Delta t}{\Delta x} \right)^2\left(u_{i+1, j}^0 -2u_{i, j}^0 + u_{i-1, j}^0\right) + \frac{1}{2}\left( \frac{c \Delta t}{\Delta y} \right)^2\left(u_{i, j+1}^0 -2u_{i, j}^0 + u_{i, j-1}^0\right).$$

In summary we have
$$
\begin{cases}
u^0_{i,j} = I(x_i, y_j), \\
u_{i,j}^{1} = u_{i, j}^0  + \Delta t V(x_i, y_j) + \frac{1}{2}\left( \frac{c \Delta t}{\Delta x} \right)^2\left(u_{i+1, j}^0 -2u_{i, j}^0 + u_{i-1, j}^0\right) + \frac{1}{2}\left( \frac{c \Delta t}{\Delta y} \right)^2\left(u_{i, j+1}^0 -2u_{i, j}^0 + u_{i, j-1}^0\right), \\
u_{i,j}^{n+1} = 2u_{i,j}^n - u_{i,j}^{n-1} + \left( \frac{c \Delta t}{\Delta x} \right)^2 (u_{i+1,j}^n - 2u_{i,j}^n + u_{i-1,j}^n) + \left( \frac{c \Delta t}{\Delta y} \right)^2 (u_{i,j+1}^n - 2u_{i,j}^n + u_{i,j-1}^n). 
\end{cases}
$$

### Boundary Conditions

We apply homogeneous Neumann boundary conditions, which specify that the normal derivative of $ u $ on the boundary is zero

$$ \frac{\partial u}{\partial n} = 0 \quad\text{for } (x, y) \in \partial \Omega, \, t \in (0, T], $$

$$
\begin{cases}
u_x(A,y,t)=\frac{u_{1, j}^n - u_{-1, j}^{n}}{2\Delta x}=0 \\
u_x(B,y,t)=\frac{u_{N_x+1, j}^n - u_{N_x-1, j}^{n}}{2\Delta x}=0 \\
u_y(x,C,t)=\frac{u_{i, 1}^n - u_{i, -1}^{n}}{2\Delta y}=0 \\
u_y(x,D,t)=\frac{u_{i, N_y+1}^n - u_{i, N_y-1}^{n}}{2\Delta y}=0 
\end{cases} \qquad \Longrightarrow \qquad

\begin{cases}
u_{1, j}^n = u_{-1, j}^{n}\\
u_{N_x+1, j}^n = u_{N_x-1, j}^{n}\\
u_{i, 1}^n = u_{i, -1}^{n}\\
u_{i, N_y+1}^n = u_{i, N_y-1}^{n}
\end{cases}.
$$

Again, to eleminate the 'fictional' values $  u_{-1, j}^{n}$, $u_{N_x+1, j}^n$, $u_{i, -1}^{n}$, $u_{i, N_y+1}^n $ we use $ (**) $. For example for $ i=0 $ we get 

$$ u_{0,j}^{n+1} = 2u_{0,j}^n - u_{0,j}^{n-1} + \left( \frac{c \Delta t}{\Delta x} \right)^2 (u_{1,j}^n - 2u_{0,j}^n + u_{-1,j}^n) + \left( \frac{c \Delta t}{\Delta y} \right)^2 (u_{0,j+1}^n - 2u_{0,j}^n + u_{0,j-1}^n), $$

$$ u_{0,j}^{n+1} = 2u_{0,j}^n - u_{0,j}^{n-1} + \left( \frac{c \Delta t}{\Delta x} \right)^2 (2u_{1,j}^n - 2u_{0,j}^n) + \left( \frac{c \Delta t}{\Delta y} \right)^2 (u_{0,j+1}^n - 2u_{0,j}^n + u_{0,j-1}^n), $$


$$ u_{0,0}^{n+1} = 2u_{0,0}^n - u_{0,0}^{n-1} + \left( \frac{c \Delta t}{\Delta x} \right)^2 (2u_{1,0}^n - 2u_{0,0}^n) + \left( \frac{c \Delta t}{\Delta y} \right)^2 (2u_{0,1}^n - 2u_{0,0}^n), $$


In the discrete setting, this is implemented by setting the boundary values to the values of their adjacent inner points:

$$ u_{0, j} = u_{1, j}, \quad u_{Nx-1, j} = u_{Nx-2, j} $$
$$ u_{i, 0} = u_{i, 1}, \quad u_{i, Ny-1} = u_{i, Ny-2} $$


