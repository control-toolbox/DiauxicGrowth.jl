# Diauxic Growth

## Context

This documentation reproduces the numerical results from Agustin G. Yabo's research on optimal control strategies for diauxic bacterial growth. The original work was implemented using [BOCOP](https://github.com/control-toolbox/bocop) (optimal control solver), and here we present a Julia implementation using the [`OptimalControl.jl`](https://control-toolbox.org/OptimalControl.jl/stable) package.

## Mathematical Framework

The mathematical model describes the time-evolution of the concentration of the ith substrate sᵢ(t), the concentration of intermediate metabolites m(t), and the volume of the cell population x(t). The states are defined as non-dimensional to simplify the computations. The dynamical system can be written as:

```math
\begin{cases}
\dot{s}_i = -w_i(s_i)\mathcal{x}, \quad i = 1, 2, \ldots, n \\
\dot{m} = \sum_{i=1}^{n} Y_i w_i(s_i) - w_R(m)(m+1) \\
\dot{\mathcal{x}} = w_R(m)\mathcal{x}
\end{cases}
```

## Optimal Control Formulation

The controlled model is:

```math
\begin{cases}
\dot{s}_i = -u_i w_i(s_i)\mathcal{x}, \quad i = 1, 2, \ldots, n \\
\dot{m} = \sum_{i=1}^{n} Y_i u_i w_i(s_i) - u_0 w_R(m)(m+1) \\
\dot{\mathcal{x}} = u_0 w_R(m)\mathcal{x}
\end{cases} \tag{S}
```

## Problem Statement

The OCP writes:

```math
\begin{cases}
\text{maximize} & J(u_0, u_1, u_2, u_3) = x(t_f) \\
\text{subject to} & \text{dynamics of (S),} \\
& \text{initial conditions: } s_i(0) = s_{i0}, m(0) = m_0, x(0) = x_0 > 0, \\
& u(\cdot) \in \mathcal{U}
\end{cases}
```

where the admissible control set is defined as:

```math
\mathcal{U} = \left\{ u = (u_0, u_1, u_2, u_3) : \begin{array}{l}
0 \leq u_i(t) \leq 1, \quad i = 0,1,2,3 \\
u_0(t) + u_1(t) + u_2(t) + u_3(t) \leq 1 \\
t \in [t_0, t_f]
\end{array} \right\}
```

## Reproducibility

```@setup main
using Pkg
using InteractiveUtils
using Markdown

# Download links for the benchmark environment
function _downloads_toml(DIR)
    link_manifest = joinpath("assets", DIR, "Manifest.toml")
    link_project = joinpath("assets", DIR, "Project.toml")
    return Markdown.parse("""
    You can download the exact environment used to build this documentation:
    - 📦 [Project.toml]($link_project) - Package dependencies
    - 📋 [Manifest.toml]($link_manifest) - Complete dependency tree with versions
    """)
end
```

```@example main
_downloads_toml(".") # hide
```

```@raw html
<details style="margin-bottom: 0.5em; margin-top: 1em;"><summary>ℹ️ Version info</summary>
```

```@example main
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📦 Package status</summary>
```

```@example main
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📚 Complete manifest</summary>
```

```@example main
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```
