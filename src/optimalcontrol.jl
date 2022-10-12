"""
Example of computing solution to optimal control problem.


References:

- https://infiniteopt.github.io/InfiniteOpt.jl/stable/examples/Optimal%20Control/consumption_savings/
- https://github.com/IlyaOrson/control_neuralode/tree/main/src
- https://discourse.julialang.org/t/solving-the-4-quadrants-of-dynamic-optimization-problems-in-julia-help-wanted/73285


"""
module OC

import InfiniteOpt, Ipopt
import Base:convert
import Distributions:Exponential, Truncated, quantile


"""
Solves the infinity horizon optimal control problem.

```{math}
\\begin{aligned}
\\max_{z,c} &  ∫_0^∞ e^{-rt} F(z(t), c(t)) dt \\\\
& s.t. \\\\
\\dot{z}(t) = & G(z(t), c(t)) \\\\
L \\leq & (c(t), z(t))' \\leq U \\\\
Z_0 = & z(0)
\\end{aligned}
```

"""
struct OCProblem{TF, TG, Tbound, Tinitial,Tr}
    F::TF
    G::TG
    L::Tbound
    U::Tbound
    dimz::Int
    dimc::Int
    z0::Tinitial
    r::Tr
end


"""
convert OCProblem to InfiniteOpt.InfininteModel
"""
function convert(::Type{InfiniteOpt.InfiniteModel}, oc::OCProblem)
    r = oc.r
    dz = oc.dimz
    dc = oc.dimc
    z0 = oc.z0
    F = oc.F 
    G = oc.G     
    m = InfiniteOpt.InfiniteModel()
    #tmass = 0.99
    #hi = quantile(Exponential(1/r),tmass)
    #t = InfiniteOpt.@infinite_parameter(m, t ~ Truncated(Exponential(1/r),0.0, hi))
    t = InfiniteOpt.@infinite_parameter(m, t ~ Exponential(1/r))
    z = InfiniteOpt.@variable(m, z[i = 1:dz], InfiniteOpt.Infinite(t))
    c = InfiniteOpt.@variable(m, c[i = 1:dc], InfiniteOpt.Infinite(t))
    for i ∈ 1:dz
        if isfinite(oc.L[i+dc]) 
            InfiniteOpt.set_lower_bound(z[i], oc.L[i+dc])
        end
        if isfinite(oc.U[i+dc]) 
            InfiniteOpt.set_upper_bound(z[i], oc.U[i+dc])
        end
    end
    for i ∈ 1:dc
        if isfinite(oc.L[i]) 
            InfiniteOpt.set_lower_bound(c[i], oc.L[i])
        end
        if isfinite(oc.U[i]) 
            InfiniteOpt.set_upper_bound(c[i], oc.U[i])
        end
    end

    
    obj = InfiniteOpt.@objective(m, Max, InfiniteOpt.@expect(F(z,c)/r,t))
    
    con = InfiniteOpt.@constraint(m,[InfiniteOpt.deriv(z[i],t) for i in 1:dz] .== G(z,c))

    for (i,z0i) ∈ enumerate(z0)
        if !isnothing(z0i)
            InfiniteOpt.@constraint(m, z[i](0) == z0i)
        end
    end
    return(m)
end

"""
    solve_InfiniteOpt(p::OCProblem; optimizer=Ipopt.Optimizer, nsupport=16)

Use InfiniteOpt.jl to solve the optimal control problem OCProblem. 

InfiniteOpt discretizes the problem on a grid of `nsupport` points. It then uses
JuMP to represent the discretized problem, and any JuMP compatible solver to maximize.

Returns a named tuple with the optimized objective value, `opt`, 
the grid of time values, `t`, grid of optimal z, `Z`, and the optimizer `status.`

Note: this does not appear to work well for the infinite horizon problems.
"""
function solve_InfiniteOpt(p::OCProblem; optimizer=Ipopt.Optimizer, nsupport=nothing)
    m = convert(InfiniteOpt.InfiniteModel, p)
    if !isnothing(nsupport)
        InfiniteOpt.fill_in_supports!(InfiniteOpt.parameter_by_name(m,"t"), num_supports = nsupport)
    end
    InfiniteOpt.set_optimizer(m, optimizer)
    InfiniteOpt.optimize!(m)
    status = InfiniteOpt.termination_status(m)    
    ts = InfiniteOpt.supports(InfiniteOpt.parameter_by_name(m,"t"))
    Z = zeros(p.dimz, length(ts))
    for i ∈ 1:p.dimz
        Z[i,:] = InfiniteOpt.value(InfiniteOpt.variable_by_name(m,"z[$i]"))
    end
    C = zeros(p.dimc, length(ts))
    for i ∈ 1:p.dimc
        C[i,:] = InfiniteOpt.value(InfiniteOpt.variable_by_name(m,"c[$i]"))
    end
    
    opt_obj = InfiniteOpt.objective_value(m)
    return(obj=opt_obj, Z=Z, C=C, t=ts, status=status)
end


import DiffEqFlux, FastGaussQuadrature, Lux, Zygote, SciMLSensitivity, Random
import DifferentialEquations: ODEFunction, ODEProblem # change to scimlbase?
import DifferentialEquations

function createode(p::OCProblem, θ, m, st;  tmax=5/p.r)
    _, re  = Lux.destructure(θ)
    function f!(dz, z, θ::AbstractArray, t)
        c = first(m([t], re(θ), st))
        dz .= p.G(z,c)        
    end
    ode_f = ODEFunction(f!)
    odeprob(a::AbstractArray) = ODEProblem(ode_f, p.z0, (0,tmax), a)    
    odeprob(tp::NamedTuple) = ODEProblem(ode_f, p.z0, (0,tmax), Lux.destructure(tp)[1])        
    return(odeprob)
end

function controlmodel(p::OCProblem; width=32)
    L = p.L
    U = p.U
    out = Vector{Function}(undef, p.dimc)
    for i ∈ 1:p.dimc
        if isfinite(L[i]) && isfinite(U[i]) 
            out[i] = x->(L[i] + (U[i]-L[i])/(1+exp(-x)))
        elseif isfinite(L[i])
            out[i] = x->(L[i] + exp(x))
        elseif isfinite(U[i])
            out[i] = x->(-exp(x) + U[i])
        else    
            out[i] = x->x
        end
    end    
    constraindomain = c->[f(x) for (f,x) ∈ zip(out,c)]

    
    m = Lux.Chain( 
            Lux.Dense(1, width, DiffEqFlux.relu),
            Lux.Dense(width, p.dimc) ,
            constraindomain)        
    return(m)
end

function solveode(odeprob; nint=100, r=0.05)
    t , ert = FastGaussQuadrature.gausslaguerre(nint)
    ert .*= r
    t .*= r    
    sol=DiffEqFlux.solve(odeprob, saveat=t)
    return(sol)
end

relaxedlogbarrier(z,δ) = (z<δ) ? (exp(1-z/δ) - 1 - log(δ)) : -log(z)


function traincontrol(θ, st, m, ocp; nint=100, iterations=100, δ=0.1,μ=1.0,
        opt= Lux.Optimisers.setup(Lux.Optimisers.ADAM(0.1f0), θ))    
    t , ert = FastGaussQuadrature.gausslaguerre(nint)
    ert .*= ocp.r
    t .*= ocp.r 
    ode = createode(ocp,θ, m, st)
    il = [i for i ∈ 1:ocp.dimz if isfinite(ocp.L[ocp.dimc+i])]
    iu = [i for i ∈ 1:ocp.dimz if isfinite(ocp.U[ocp.dimc+i])]
    penalty(z) = μ*(sum(relaxedlogbarrier.(z[i] - ocp.L[ocp.dimc+i], δ) for i ∈ il) + 
        sum(relaxedlogbarrier.(-z[i] + ocp.U[ocp.dimc+i], δ) for i ∈ iu)) 
    function loss(m, θ, st) 
        sol = DiffEqFlux.solve(ode(Lux.destructure(θ)[1]),saveat=t
            #,sensealg=SciMLSensitivity.QuadratureAdjoint()
            )[1,:]
        l = sum(-ocp.F(s, first(m([τ], θ, st)))*w + penalty(s) for (τ, w, s) ∈ zip(t,ert,sol))
        return(l)
    end
      
    ℓ = p -> loss(m, p, st)
    lval, back = Zygote.pullback(ℓ, θ)
    gs = back(one(lval))[1]
    @show gs
    for i ∈ 1:iterations
        lval, back = Zygote.pullback(ℓ, θ)
        gs = back(one(lval))[1]
        opt, θ = Lux.Optimisers.update(opt, θ, gs)

        println("Iter[$i]: obj=$(-lval)")
    end
    return(θ, m, st, opt)
end


end