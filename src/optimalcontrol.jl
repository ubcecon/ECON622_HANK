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
import Distributions:Exponential


"""
Solves the infinity horizon optimal control problem.

```{math}
\\begin{aligned}
\\max_z &  ∫_0^∞ e^{-rt} F(z(t),t) dt \\\\
& s.t. \\\\
0 = & G(z(t), \\dot{z}(t), t) \\\\
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

function convert(::Type{InfiniteOpt.InfiniteModel}, oc::OCProblem)
    r = oc.r
    d = oc.dim
    z0 = oc.z0
    F = oc.F 
    G = oc.G     
    m = InfiniteOpt.InfiniteModel()
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


function solve_InfiniteOpt(p::OCProblem; optimizer=Ipopt.Optimizer, nsupport=16)
    m = convert(InfiniteOpt.InfiniteModel, p)
    InfiniteOpt.fill_in_supports!(InfiniteOpt.parameter_by_name(m,"t"), num_supports = nsupport);
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

end
