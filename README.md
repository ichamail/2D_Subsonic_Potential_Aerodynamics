# Subsonic Potential Aerodynamics

## A 2D Aerodynamic Potential-Flow Code

### Theoretical Model

![Potential Flow Theoretical Model](https://github.com/ichamail/2D-Subsonic-Potential-Aerodynamics/assets/107580530/8ddf6ab4-904d-4ab3-a8a4-5274766b0f3b)


#### Velocity Field
   * incompressible: $` \nabla \cdot \underline{V} = 0 `$
   * conservative: $` \underline{V} = \nabla \phi \implies \nabla \times \underline{V} = \nabla \times \nabla \phi = 0 `$ (irrotational)

```math
\nabla \cdot \underline{V} = \nabla \cdot \nabla \phi = \nabla^2 \phi = 0
```


#### Integral Equation of velocity potential $\phi$
```math
\phi(\underline{r}_p) =
\iint_S \frac{\sigma}{2\pi} \ln{(\lVert \underline{r} - \underline{r}_p \rVert)} \mathrm{d}S + \iint_{S \cup S_w}  - \frac{\mu}{2\pi} (\underline{e}_n \cdot \nabla) \ln{(\lVert \underline{r} - \underline{r}_p \rVert)} \mathrm{d}S + \phi_\infty(\underline{r}_p)
```


where:
   * $` \mu = \phi - \phi_i `$
   * $` \sigma = (\underline{e}_n \cdot \nabla)(\phi - \phi_i) `$
   * $` \phi_\infty(\underline{r}_p) = \iint_{S_\infty}  \left[ (\underline{e}_n \cdot \nabla) \phi \ln{ \left( \frac{\lVert \underline{r} - \underline{r}_p \rVert}{2\pi} \right) } - \phi  (\underline{e}_n \cdot \nabla) \ln{ \left( \frac{\lVert \underline{r} - \underline{r}_p \rVert}{2\pi} \right) } \right] \mathrm{d}S `$


if $` \mu = \phi - \phi_i = \phi - \phi_\infty `$ and $` \sigma = (\underline{e}_n \cdot \nabla)(\phi - \phi_i) = (\underline{e}_n \cdot \nabla)(\phi - \phi_\infty ) `$, then for $` P \in S^- `$ :

```math
\iint_S \frac{\sigma}{2\pi} \ln{(\lVert \underline{r} - \underline{r}_p \rVert)} \mathrm{d}S + \iint_{S \cup S_w} - \frac{\mu}{2\pi} (\underline{e}_n \cdot \nabla) \ln{(\lVert \underline{r} - \underline{r}_p \rVert)} \mathrm{d}S = 0, \qquad \forall (x_p, y_p, z_p) \in S: (\underline{r} - \underline{r}_p) \cdot \underline{e}_n \to 0^+
```

note: $` \sigma = (\underline{e}_n \cdot \nabla)(\phi - \phi_i) = (\underline{e}_n \cdot \nabla)(\phi - \phi_\infty ) = (\underline{e}_n \cdot \nabla)\phi - (\underline{e}_n \cdot \nabla)\phi_\infty = \underline{V} - \underline{e}_n \cdot \underline{V}_\infty = - \underline{e}_n \cdot \underline{V}_\infty `$

Notable Observations:
   * A vector field $` \underline{V}: U \to \mathbb{R}^n `$, where $` \mathbb{U} `$ is an open subset of $` \mathbb{R}^n `$, is said to be conservative if  there exists a $` \mathrm{C}^1 `$ (continuously differentiable) scalar field $` \phi `$ on $` \mathbb{U} `$ such that $` \underline{V} = \nabla \phi `$.

   * According to Poincaré's Lemma, A continuously differentiable ($` \mathrm{C}^1 `$) vector field $` \underline{V} `$ defined on a simply connected subset $` \mathbb{U} `$ of $` \mathbb{R}^n `$  ($` \underline{V} \colon \mathbb{U} \subseteq \mathbb{R}^n \to \mathbb{R}^n `$), is conservative if and only if it is irrotational throughout its domain ($` \nabla \times \underline{V} = 0 `$, $` \forall \underline{x} \in \mathbb{U} `$).

   * Circulation $` \Gamma = \oint_{C} \underline{V} \cdot \mathrm{d} \underline{l} = \iint_S \nabla \times \underline{V} \cdot \mathrm{d}\underline{S} `$ 
   In a conservative vector field this integral evaluates to zero for every closed curve. $` \Gamma = \oint_{C} \underline{V} \cdot \mathrm{d} \underline{l} = \iint_S \nabla \times \underline{V} \cdot \mathrm{d}\underline{S} = \iint_S \nabla \times \nabla \phi \cdot \mathrm{d}\underline{S} = 0 `$

   * The space exterior to a 2D object (e.g. an airfoil) is not simply connected

### Numerical Model (Panel Methods)

```math
\sum_{j=0}^{N_s - 1} B_{ij} \sigma_j + \sum_{j=0}^{N_s + N_w - 1} C_{ij} \mu_j = 0 , \qquad 0 \le i < N_s
``` 

where:
   * $` B_{ij} =  \frac{1}{2\pi} \iint_{S_j}  \ln{(\lVert \underline{r} - \underline{r}_{cp_i} \rVert)} \mathrm{d}{S_j}  = \frac{1}{2\pi} \iint_{t_1}^{t_2}  \ln \left(\sqrt{(t_{cp_i} - t)^2 + n_{cp_i}^2} \right) \mathrm{d}{t_j} `$

   * $` C_{ij} =  \frac{1}{2\pi} \iint_{S_j}  (\underline{e}_n \cdot \nabla) \ln{(\lVert \underline{r} - \underline{r}_{cp_i} \rVert)} \mathrm{d}{S_j} = \frac{1}{2\pi} \iint_{S_j}  \frac{n_{cp_i}}{(t_{cp_i} - t)^2 + n_{cp_i}^2} \mathrm{d}{t_j} `$

from Kutta Condition: $` \mu_w = const = \mu_U - \mu_L `$
```math
A_{ij} \mu_j = - B_{ij} \sigma_j , \qquad A_{ij} = 
\begin{cases}
   C_{ij} + \sum\limits_{k} C_{ik} & \text{if $j=0$}\\
   C_{ij} & \text{if $0 < j < N_s - 1$}\\
   C_{ij} - \sum\limits_{k} C_{ik} & \text{if $j=N_s-1$}
\end{cases} 
```


```math
0 \le i < N_s  \qquad 0 \le j < N_s  \qquad N_s \le k < N_s + N_w
```

### Features
 1. Calculation of Non-Lifting Potential Flow about 2D arbitrarily-shaped rigid bodies
    1. Steady simulations
    2. Unsteady simulations

 2. Calculation of Lifting Pseudo-Potential Flow around 2D arbitrarily-shaped rigid bodies
      1. Steady state simulations with flat rigid wake model
      2. Steady state iterative simulations with flexible wake model 
      3. Unsteady simulations with a shedding wake model

## Simulation Results

## Potential Flow around a 2D Circular Object


![circle mesh in body-fixed frame](/images/potential_flow_around_2D_circular_object/circle_mesh_in_inertial_frame_of_reference.png)

![circle mesh in body-fixed frame](/images/potential_flow_around_2D_circular_object/circle_mesh_in_body-fixed_frame_of_reference.png)

![circle stream plot 10 panels](/images/potential_flow_around_2D_circular_object/circle_streamplot_10_panels.png)

![circle pressure coefficient 10 panel](/images/potential_flow_around_2D_circular_object/circle_pressure_coefficient_10_panels.png)

![circle pressure coefficient 20 panel](/images/potential_flow_around_2D_circular_object/circle_pressure_coefficient_20_panels.png)


## Potential Flow around an Airfoil

![numerical model in inertial frame](/images/potential_flow_around_airfoil/steady_simulation/numerical_model_inertial_frame.png)

![numerical model in body-fixed frame](/images/potential_flow_around_airfoil/steady_simulation/numerical_model_bodyfixed_frame.png)

## Steady Simulation with rigid wake

![rigid wake](/images/potential_flow_around_airfoil/steady_simulation/rigid_wake.png)



![stream lines](/images/potential_flow_around_airfoil/steady_simulation/streamlines.png)

![pressure coefficient distribution](/images/potential_flow_around_airfoil/steady_simulation/pressure_coefficient_distribution.png)

![aerodynamic force vector](/images/potential_flow_around_airfoil/steady_simulation/aerodynamic_force_vector.png)



## Steady Simulation with iterative wake

![iterative wake](/images/potential_flow_around_airfoil/steady_iterative_simulation/iterative_wake.png)

![streamlines](/images/potential_flow_around_airfoil/steady_iterative_simulation/streamlines.png)


![pressure coefficient distribution](/images/potential_flow_around_airfoil/steady_iterative_simulation/pressure_coefficient_distribution.png)

![aerodynamic force vector](/images/potential_flow_around_airfoil/steady_iterative_simulation/aerodynamic_force_vector.png)




## Unsteady Simulation with wake roll up

![wake rollup](/images/potential_flow_around_airfoil/unsteady_simulation/wake_rollup.png)

![streamlines](/images/potential_flow_around_airfoil/unsteady_simulation/streamlines.png)

![wake rollup final step](/images/potential_flow_around_airfoil/unsteady_simulation/wake_rollup_final_step.png)

![pressure coefficient distribution](/images/potential_flow_around_airfoil/unsteady_simulation/pressure_coefficient_distribution.png)

![aerodynamic force vector](/images/potential_flow_around_airfoil/unsteady_simulation/aerodynamic_force_vector.png)



