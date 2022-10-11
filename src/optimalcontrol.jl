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
\\max_z &  ∫_0^∞ e^{-rt} F(z(t)) dt \\\\
& s.t. \\\\
0 = & G(z(t), \\dot{z}(t)) \\\\
L \\leq & z(t) \\leq U \\\\
Z_0 = & z(0)
\\end{aligned}
```

"""
struct OCProblem{TF, TG, Tbound, Tinitial,Tr}
    F::TF
    G::TG
    L::Tbound
    U::Tbound
    dim::Int
    z0::Tinitial
    r::Tr
end


"""
convert OCProblem to InfiniteOpt.InfininteModel
"""
function convert(::Type{InfiniteOpt.InfiniteModel}, oc::OCProblem)
    r = oc.r
    d = oc.dim
    z0 = oc.z0
    F = oc.F 
    G = oc.G     
    m = InfiniteOpt.InfiniteModel()
    #tmass = 0.99
    #hi = quantile(Exponential(1/r),tmass)
    #t = InfiniteOpt.@infinite_parameter(m, t ~ Truncated(Exponential(1/r),0.0, hi))
    t = InfiniteOpt.@infinite_parameter(m, t ~ Exponential(1/r))
    z = InfiniteOpt.@variable(m, z[i = 1:d], InfiniteOpt.Infinite(t))
    for i ∈ 1:d
        if isfinite(oc.L[i]) 
            InfiniteOpt.set_lower_bound(z[i], oc.L[i])
        end
        if isfinite(oc.U[i]) 
            InfiniteOpt.set_upper_bound(z[i], oc.U[i])
        end
    end

    
    obj = InfiniteOpt.@objective(m, Max, InfiniteOpt.@expect(F(z)/r,t))
    
    con = InfiniteOpt.@constraint(m, G(z,[InfiniteOpt.deriv(z[i],t) for i ∈ 1:d])==0)

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
    Z = zeros(p.dim, length(ts))
    for i ∈ 1:p.dim
        Z[i,:] = InfiniteOpt.value(InfiniteOpt.variable_by_name(m,"z[$i]"))
    end
    opt_obj = InfiniteOpt.objective_value(m)
    return(obj=opt_obj, Z=Z, t=ts, status=status)
end


import DiffEqFlux, FastGaussQuadrature, Lux, Zygote, SciMLSensitivity, Random
import DifferentialEquations: DAEFunction, DAEProblem # change to scimlbase?
import DifferentialEquations

function createode(p::OCProblem, m, st; controls=1:(p.dim÷2), states=setdiff(1:p.dim, controls), tmax=5/p.r)
    θ = Lux.initialparameters(Random.default_rng(), m)
    _, re  = Lux.destructure(θ)
    function f!(resid, du, u, θ, t)
        z = Vector{eltype(u)}(undef,p.dim)
        dotz = Vector{eltype(du)}(undef,p.dim)
        z[controls] .= first(m([t], re(θ), st))
        z[states] .= u
        dotz[controls] .= NaN
        dotz[states] .= du
        resid .= p.G(z,dotz)        
    end
    dae_f = DAEFunction(f!)
    odeprob(θ::AbstractArray) = DAEProblem(dae_f, [0.0], p.z0[states], (0,tmax), p=θ)    
    odeprob(θ::NamedTuple) = DAEProblem(dae_f, [0.0], p.z0[states], (0,tmax), p=Lux.destructure(θ)[1])        
    return(odeprob)
end

function controlmodel(p::OCProblem; controls=1:(p.dim ÷ 2), width=32)
    L = p.L[controls]
    U = p.U[controls]
    function maketransform(c)
        out = Vector{Function}(undef, length(c))
        for i ∈ eachindex(c)
            if isfinite(L[i]) && isfinite(U[i]) 
                out[i] = x->(L[i] + (U[i]-L[i])/(1+exp(x)))
            elseif isfinite(L[i])
                out[i] = x->(L[i] + exp(x))
            elseif isfinite(U[i])
                out[i] = x->(-exp(x) + U[i])
            else    
                out[i] = x->x
            end
        end

        return(c->[f(x) for (f,x) ∈ zip(out,c)])
    end
    constraindomain = maketransform(controls)

    
    m = Lux.Chain( 
            Lux.Dense(1, width, DiffEqFlux.relu),
            Lux.Dense(width, length(controls)) ,
            constraindomain)        
    return(m)
end

function solveode(odeprob, θ; nint=100, r=0.05)
    t , ert = FastGaussQuadrature.gausslaguerre(nint)
    ert .*= r
    t .*= r    
    sol=DifferentialEquations.solve(odeprob, DifferentialEquations.DImplicitEuler(), p=θ,saveat=t)
    return(sol)
end
solveode(p::OCProblem,θ, m, st) = solveode(createode(p, m, st)(θ))

import Base: length
struct Control{T}
    m::T
end
(c::Control)(t) = c.m([t])
length(c::Control) = 100+sum(length(s) for s ∈ Flux.params(c.m)) # any number > 100 to trigger reverse mode 

function traincontrol!(θ, st, m, ocp; nint=100, iterations=100)    
    t , ert = FastGaussQuadrature.gausslaguerre(nint)
    ert .*= ocp.r
    t .*= ocp.r 
    odesolver = DifferentialEquations.DImplicitEuler()
    ode = createode(ocp,m, st)(θ)
    function loss(m, θ, st) 
        sol = DifferentialEquations.solve(ode, odesolver, p=Lux.destructure(θ)[1],saveat=t
            #,sensealg=SciMLSensitivity.QuadratureAdjoint()
            )[1,:]
        #sum(ocp.F(vcat(m([τ]), sol(τ)))*w for (τ,w) ∈ zip(t, ert))        
        #sum(sum(vcat(m([τ])))*w for (τ,w,s) ∈ zip(t, ert, sol))
        l = sum(-ocp.F(vcat(first(m([s], θ, st)),s))*w for (τ, w, s) ∈ zip(t,ert,sol))
        return(l)
    end
    opt= Lux.Optimisers.setup(Lux.Optimisers.ADAM(0.1f0), θ)
    #opt = DiffEqFlux.ADAM(0.1)
    #cb = function ()
    #    display(loss())
    #end
    
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
    return(θ, m, st)
end


end