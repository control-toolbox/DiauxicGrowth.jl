# Diauxic Growth

## Context

This documentation reproduces the numerical results from Agustin G. Yabo's research on optimal control strategies for diauxic bacterial growth. The original work was implemented using BOCOP (optimal control solver), and here we present a Julia implementation using the `OptimalControl.jl` package.

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

### Problem Statement

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

## Implementation and Results

```@example diauxic
using DifferentialEquations
using Plots
using OptimalControl
using NLPModelsIpopt
using LaTeXStrings
using Plots.PlotMeasures

```

```@example diauxic
# Define parameters first
k = 1
KR = 0.003
Kᵢ = 0.0003

# Define uptake functions
wm1(s) = k*s/(Kᵢ + s)
wm2(s) = k*s/(Kᵢ + s)
wm3(s) = k*s/(Kᵢ + s)
wR(m) = k * m / (KR + m)

function diauxci_ocp(t₀, tf, I₀, Y₁, Y₂, Y₃)
    ocp = @def begin
        # time interval
        t ∈ [t₀, tf], time
        # state variables: z = (s₁, s₂, s₃, m, x)  
        z = (s₁, s₂, s₃, m, x) ∈ R⁵, state   
        # control variables: u = (u₀, u₁, u₂, u₃)
        u = (u₀, u₁, u₂, u₃) ∈ R⁴, control   
        # Control bounds: each u_i ∈ [0,1]
        [0, 0, 0, 0] ≤ u(t) ≤ [1, 1, 1, 1]    
        u₀(t) + u₁(t) + u₂(t) + u₃(t) ≤ 1 # Control bounds: sum of controls ≤ 1

        z(t) ≥ [0, 0, 0, 0, 0]           # State bounds: all states non-negative

        z(t₀) == [I₀[1], I₀[2], I₀[3], I₀[4], I₀[5]]  # Initial conditions

        ż(t) == [
            - u₁(t) * wm1(s₁(t)) * x(t),                          
            - u₂(t) * wm2(s₂(t)) * x(t),                            
            - u₃(t) * wm3(s₃(t)) * x(t),                          
              Y₁*u₁(t)*wm1(s₁(t)) + Y₂*u₂(t)*wm2(s₂(t))         
              + Y₃*u₃(t)*wm3(s₃(t)) - u₀(t)*wR(m(t))*(m(t)+1),                                 
              u₀(t)*wR(m(t)) * x(t)            
        ]
        x(tf) → max     # Maximize final biomass
    end 
    
    sol = solve(
        ocp, :direct, :adnlp, :ipopt;
        disc_method = :gauss_legendre_3,
        grid_size   = 1000,       
        tol         = 1e-8,
        display     = true)

    return sol 
end
```

```@example diauxic
# Initial conditions and time horizon
s₁₀, s₂₀, s₃₀, m₀, x₀ = 1e-3, 2e-3, 3e-3, 1e-3, 5e-3 
t₀, tf = 0.0, 2.0
Y₁, Y₂, Y₃ = 1.0, 0.5, 0.2
I₀ = [s₁₀, s₂₀, s₃₀, m₀, x₀]

# Solve the optimal control problem
sol = diauxci_ocp(t₀, tf, I₀, Y₁, Y₂, Y₃)
```

## Switching Time Analysis

```@example diauxic
function calculate_switching_times(sol, Y₁, Y₂, Y₃)
    """
    Calculate switching times based on when yield-weighted Michaelis-Menten functions
    
    Returns:
    - arc1: Time when Y₁*wm1 = Y₂*wm2 
    - arc2: Time when Y₂*wm2 = Y₃*wm3
    """
    arc1, arc2 = nothing, nothing
    
    try 
        t = time_grid(sol)
        z = state(sol)
        
        # Get state values at each time point
        s1_vals = [z(ti)[1] for ti in t]
        s2_vals = [z(ti)[2] for ti in t]
        s3_vals = [z(ti)[3] for ti in t]
        
        # Evaluate the Michaelis-Menten functions
        wm1_vals = [wm1(s1) for s1 in s1_vals]
        wm2_vals = [wm2(s2) for s2 in s2_vals]
        wm3_vals = [wm3(s3) for s3 in s3_vals]
        
        # Find switching times where yield-weighted functions are equal
        arc1_idx = findfirst((Y₁ .* wm1_vals .- Y₂ .* wm2_vals).^2 .< 1e-3)
        arc2_idx = findfirst((Y₂ .* wm2_vals .- Y₃ .* wm3_vals).^2 .< 1e-3)
        
        arc1 = t[arc1_idx]
        arc2 = t[arc2_idx]    
    catch e
        println("Error in switching time calculation: ", e)
    end
    
    return arc1, arc2
end

# Calculate switching times for our solution
arc1, arc2 = calculate_switching_times(sol, Y₁, Y₂, Y₃)
```

## State Variable Visualization

```@example diauxic
function plot_state_variables(sol, t₀, tf)
    """
    Generate three separate plots for the state variables:
    1. Substrates (s₁, s₂, s₃)
    2. Metabolites (m)
    3. Biomass (x)
    
    Returns a tuple of three plots: (substrates_plot, metabolites_plot, biomass_plot)
    """
    
    # Plot 1: Substrates
    substrates_plot = plot(
        t -> state(sol)(t)[3], 0, tf,
        label = L"$s_3$", color = "#79af97", lw = 3,
        grid = true,
        gridlinewidth = 1.5,
        size = (400, 400),
        legendfontsize = 10,
        tickfontsize = 6,
        guidefontsize = 10,
        framestyle = :box,
        left_margin = 10mm,
        bottom_margin = 10mm,
        foreground_color_subplot = :black, 
        xlabel = "t",
        ylabel = "Extracellular quantities"
    )
    plot!(substrates_plot,
        t -> state(sol)(t)[2], 0, tf,
        label = L"$s_2$", color = "#df8f44", lw = 3
    )
    plot!(substrates_plot,
        t -> state(sol)(t)[1], 0, tf,
        label = L"$s_1$", color = "#b24746", lw = 3
    )

    # Plot 2: Metabolites
    metabolites_plot = plot(
        t -> state(sol)(t)[4], 0, tf,
        label = L"$m$", color = "#374e55", lw = 3,
        grid = true,
        gridlinewidth = 1.5,
        size = (300, 300),
        legendfontsize = 7,
        tickfontsize = 10,
        guidefontsize = 13,
        framestyle = :box,
        left_margin = 10mm,
        bottom_margin = 10mm,
        foreground_color_subplot = :black, 
        xlabel = "t",
        ylabel = "metabolite"
    )

    # Plot 3: Biomass
    biomass_plot = plot(
        t -> state(sol)(t)[5], 0, tf,
        label = L"$x$", color = "#374e55", lw = 3,
        grid = true,
        gridlinewidth = 1.5,
        size = (300, 300),
        legendfontsize = 7,
        tickfontsize = 10,
        guidefontsize = 13,
        framestyle = :box,
        left_margin = 10mm,
        bottom_margin = 10mm,
        foreground_color_subplot = :black, 
        xlabel = "t",
        yformatter = :plain,
        ylabel = "biomass"
    )
    
    return substrates_plot, metabolites_plot, biomass_plot
end

# Generate the three state variable plots
substrates_plot, metabolites_plot, biomass_plot = plot_state_variables(sol, t₀, tf)
```
```@example diauxic
substrates_plot
```

```@example diauxic
metabolites_plot
```

```@example diauxic
biomass_plot
```

## Control Analysis

```@example diauxic
function plot_controls(sol, t₀, tf; arc1=nothing, arc2=nothing)
    """
    Plot optimal controls as stacked areas with switching time indicators
    
    Arguments:
    - sol: Solution from optimal control problem
    - t₀, tf: Time bounds
    - arc1, arc2: Optional switching times to highlight with vertical lines
    """
    
    # Extract time grid and control function
    t_vals = time_grid(sol)
    u = control(sol)
    
    # Evaluate controls at time points
    u0_vals = [u(t)[1] for t in t_vals]
    u1_vals = [u(t)[2] for t in t_vals]
    u2_vals = [u(t)[3] for t in t_vals] 
    u3_vals = [u(t)[4] for t in t_vals]
    
    # Find switching times 
    arcr1_idx = findfirst(u0_vals .> 0.1)
    arcr2_idx = findfirst(u0_vals .> 0.5)
    
    arcr1 = arcr1_idx !== nothing ? t_vals[arcr1_idx] : nothing
    arcr2 = arcr2_idx !== nothing ? t_vals[arcr2_idx] : nothing
				   
    
    p = plot(t_vals, u0_vals, 
         fillrange=0, fillalpha=0.6, fillcolor="#374e55",
         label=L"$u_0$", linewidth=0, 
         grid=true, gridlinewidth=1.5,
         size=(1600, 1000),
         legendfontsize=24, tickfontsize=20, guidefontsize=26,
         framestyle=:box, left_margin=10mm, bottom_margin=10mm,
         foreground_color_subplot=:black,
         xlabel="t", ylabel="Optimal controls",
         xlims=(t₀, tf), ylims=(0, 1)
    )
    
    # u1 on top of u0
    plot!(p, t_vals, u0_vals .+ u1_vals, 
          fillrange=u0_vals, fillalpha=0.9, fillcolor="#b24746", 
          label=L"$u_1$", linewidth=0)
    
    # u2 on top of u0+u1  
    plot!(p, t_vals, u0_vals .+ u1_vals .+ u2_vals, 
          fillrange=(u0_vals .+ u1_vals), fillalpha=0.9, fillcolor="#df8f44", 
          label=L"$u_2$", linewidth=0)
    
    # u3 on top of u0+u1+u2
    plot!(p, t_vals, u0_vals .+ u1_vals .+ u2_vals .+ u3_vals, 
          fillrange=(u0_vals .+ u1_vals .+ u2_vals), fillalpha=0.9, 
          fillcolor="#79af97", label=L"$u_3$", linewidth=0)
    
    plot!(p, t_vals, u0_vals, linewidth=0.5, linecolor=:black, alpha=1, label="")
    plot!(p, t_vals, u0_vals .+ u1_vals, linewidth=0.5, linecolor=:black, label="")
    plot!(p, t_vals, u0_vals .+ u1_vals .+ u2_vals, linewidth=0.5, linecolor=:black, label="")
    
    # Adding vertical lines for switching times
    if arc1 !== nothing
        vline!(p, [arc1], color=:black, linestyle=:dash, linewidth=1.5, alpha=0.9, label="")
    end
    if arc2 !== nothing
        vline!(p, [arc2], color=:black, linestyle=:dash, linewidth=1.5, alpha=0.9, label="")
    end
    if arcr1 !== nothing
        vline!(p, [arcr1], color=:black, linestyle=:dash, linewidth=1.5, alpha=0.9, label="")
    end
    if arcr2 !== nothing
        vline!(p, [arcr2], color=:black, linestyle=:dash, linewidth=1.5, alpha=0.9, label="")
    end

    plot!(p, legend=:topright, legend_columns=2)
    
    return p
end

# Plot the controls for our solution
control_plot = plot_controls(sol, t₀, tf; arc1=arc1, arc2=arc2)
```

## References

The mathematical models and optimal control strategies implemented in this documentation are based on the following research:

1. **Yabo, A. G.** (2023). Predicting microbial cell composition and diauxic growth as optimal control strategies. *IFAC-PapersOnLine*, 56(2), 6217-6222. [https://doi.org/10.1016/j.ifacol.2023.10.745.](https://www.sciencedirect.com/science/article/pii/S2405896323011229)

2. **Yabo, A. G.** (2025). Optimal control strategies in a generic class of bacterial growth models with multiple substrates. *Automatica*, 171, 111881. [https://doi.org/10.1016/j.automatica.2024.111881](https://doi.org/10.1016/j.automatica.2024.111881)

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
