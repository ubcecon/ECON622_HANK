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


hamiltonian(p) = (t,z,c,λ)->exp(-p.r*t)*p.F(z,c) + λ'*p.G(z,c)


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
    tmass = 0.99
    hi = quantile(Exponential(1/r),tmass)
    #t = InfiniteOpt.@infinite_parameter(m, t ~ Truncated(Exponential(1/r),0.0, hi))
    t = InfiniteOpt.@infinite_parameter(m, t ∈ [0, hi])
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

    
    obj = InfiniteOpt.@objective(m, Max, InfiniteOpt.@integral(F(z,c)*exp(-r*t),t))
    
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


import FastGaussQuadrature, Lux, Zygote, SciMLSensitivity, Random
import CommonSolve
import SciMLBase: ODEFunction, ODEProblem 
import OrdinaryDiffEq: Tsit5


evalm(t,z::AbstractArray,m,θ,st) = first(m(z, θ, st))
evalm(t,z,m,θ,st) = first(m([z], θ, st))

function createode(p::OCProblem, θ, m, st;  tmax=5/p.r)
    _, re  = Lux.destructure(θ)
    function f!(dz, z, θ::AbstractArray, t)
        c = evalm(t,z,m, re(θ),st)
        dz .= p.G(z,c)        
    end
    ode_f = ODEFunction(f!)
    odeprob(a::AbstractArray) = ODEProblem(ode_f, p.z0, (0,tmax), a)    
    odeprob(tp::NamedTuple) = ODEProblem(ode_f, p.z0, (0,tmax), Lux.destructure(tp)[1])        
    return(odeprob)
end

function constraindomain(L,U, T=eltype(L))
    out = Vector{Function}(undef, length(L))
    L = T.(L)
    U = T.(U)
    for i ∈ 1:length(L)
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
    return(c->[f(x) for (f,x) ∈ zip(out,c)])
end

function controlmodel(p::OCProblem; width=32, maxt=80.0/p.r)
    @views cons = constraindomain(p.L[1:p.dimc], p.U[1:p.dimc])
    scalet = t->t/(maxt)*2
    m = Lux.Chain( 
            scalet,
            Lux.Dense(1, width, Lux.relu),
            Lux.Dense(width, p.dimc) ,
            cons)        
    return(m)
end

function solveode(odeprob; nint=100, r=0.05)
    t , ert = FastGaussQuadrature.gausslaguerre(nint)
    ert ./= r
    t ./= r    
    sol=CommonSolve.solve(odeprob, Tsit5(), saveat=t)
    return(sol)
end

relaxedlogbarrier(z,δ) = (z<δ) ? (exp(1-z/δ) - 1 - log(δ)) : -log(z)


function traincontrol(θ, st, m, ocp; nint=100, iterations=100, δ=0.1,μ=1.0,
        opt= Lux.Optimisers.setup(Lux.Optimisers.ADAM(0.1f0), θ))    
    t , ert = FastGaussQuadrature.gausslaguerre(nint)
    ert ./= ocp.r
    t ./= ocp.r 
    ode = createode(ocp,θ, m, st)
    il = [i for i ∈ 1:ocp.dimz if isfinite(ocp.L[ocp.dimc+i])]
    iu = [i for i ∈ 1:ocp.dimz if isfinite(ocp.U[ocp.dimc+i])]
    penalty(z) = μ*(sum(relaxedlogbarrier.(z[i] - ocp.L[ocp.dimc+i], δ) for i ∈ il) + 
        sum(relaxedlogbarrier.(-z[i] + ocp.U[ocp.dimc+i], δ) for i ∈ iu)) 
    function loss(m, θ, st) 
        sol = CommonSolve.solve(ode(Lux.destructure(θ)[1]),Tsit5(), saveat=t, maxiters=1e3
            #,sensealg=SciMLSensitivity.QuadratureAdjoint()
            )[1,:]
        l = sum(-ocp.F(z, evalm(τ, z, m, θ, st))*w + penalty(z) for (τ, w, z) ∈ zip(t,ert,sol))
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

import FiniteDiff
"""

Method of Mortezaee & Nazemi (2019)
"""
function neuralcollocation(p, width; points=100, opt = nothing)
    RType = Float32
    t , ert = FastGaussQuadrature.gausslaguerre(points)
    ert ./= p.r
    t ./= p.r

    invmaxt = 1/maximum(t)
    scalet = t->t*invmaxt*2
    cm = Lux.Chain( scalet,
             Lux.Dense(1, width, Lux.relu),
            Lux.Dense(width, p.dimc) ,
        constraindomain(p.L[1:p.dimc], p.U[1:p.dimc], RType))        
    zm = Lux.Chain( scalet,
        Lux.Dense(1, width, Lux.relu),
        Lux.Dense(width, p.dimz) ,
        constraindomain(p.L[p.dimc+1:end].-p.z0, p.U[p.dimc+1:end].-p.z0, RType))        
    λm = Lux.Chain( scalet,
            Lux.Dense(1, width, Lux.relu),
            Lux.Dense(width, p.dimz))                    

    rng = Random.default_rng()
    θc, stc = Lux.setup(rng,cm)
    θc.layer_2.bias .= Random.rand(scalet.(t), width)
    θc.layer_3.weight .*= 0.1
    θc.layer_3.bias .= -log((p.U[1]-p.L[1])/(1.0+p.r*p.z0[1]-p.L[1])-1) 
    θz, stz = Lux.setup(rng,zm)
    θz.layer_2.bias .= Random.rand(scalet.(t), width)
    θz.layer_3.bias .= -log((p.U[2]-p.L[2])/(p.z0[1]-p.L[2])-1) 
    θλ, stλ = Lux.setup(rng,λm)
    θλ.layer_2.bias .= Random.rand(scalet.(t), width)

    H = hamiltonian(p)
    focs(t, z, c, λ, ż, λ̇) = [ 
        ż .- p.G(z,c),
        FiniteDiff.finite_difference_gradient(c->H(t,z,c,λ),c,Val(:central), eltype(c), Val(false)),
        #ForwardDiff.gradient(c->H(t,z,c,λ),c),
        #λ̇ + FiniteDiff.finite_difference_gradient(z->H(t,z,c,λ),z,Val(:central), eltype(z), Val(false))
           # ForwardDiff.gradient(z->H(t,z,c,λ),z)
        ]
    RType = eltype(θc.layer_2.bias)
    t = RType.(t)
    function loss(θc, θz, θλ)
        C(t) = cm([t], θc, stc)[1]
        Z(t) = RType.(p.z0) .+ tanh(t)*zm([t], θz, stz)[1]
        Ż(t) = FiniteDiff.finite_difference_derivative(Z,t)
        Λ(t) = exp(-RType(p.r)*t)*λm([t], θλ, stλ)[1]
        Λ̇(t) = FiniteDiff.finite_difference_derivative(Λ,t)
        sum( sum(x->x'*x, 
            focs(s, Z(s), C(s), Λ(s), Ż(s), Λ̇(s)))
         for s ∈ t
        )
    end
    
    
    ℓ = p -> loss(p...)
    θ = (θc,θz,θλ)
    lval, back = Zygote.pullback(ℓ,θ)
    gs = back(one(lval))[1]
    @show gs
    x, re = Lux.destructure(θ)
    gf = ForwardDiff.gradient(x->ℓ(re(x)),x)
    
    for (z, f) ∈ zip(gs, re(gf))        
        for k ∈ keys(z)
            if !(k ∈ keys(f))
                @warn "$k not found in forward gradient"
                continue
            end
            try
                @show k
                @show norm(z[k][:bias]-f[k][:bias])
                @show norm(z[k][:weight]-f[k][:weight])
            catch
                println("$k doesn't have bias and/or weight")
            end
        end
    end

    if isnothing(opt)
        opt = Lux.Optimisers.setup(Lux.Optimisers.ADAM(0.1f0), θ)
    end 
    for i ∈ 1:iterations
        lval, back = Zygote.pullback(ℓ, θ)
        gs = back(one(lval))[1]
        opt, θ = Lux.Optimisers.update(opt, θ, gs)
        println("Iter[$i]: obj=$(lval)")
    end
    θc, θz, θλ = θ
    C(t) = cm([t], θc, stc)[1]
    Z(t) = p.z0 + tanh(t)*zm([t], θz, stz)[1]
    Ż(t) = ForwardDiff.derivative(Z,t)
    Λ(t) = exp(p.r*t)*λm([t], θλ, stλ)[1]
    Λ̇(t) = ForwardDiff.derivative(Λ,t)
    return(C=C, Z=Z, Ż = Ż, Λ=Λ, Λ̇=Λ̇, θ=θ, opt=opt)    
end

end