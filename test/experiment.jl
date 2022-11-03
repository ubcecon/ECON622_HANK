import HANK, FastGaussQuadrature, PlotlyLight, Lux, Zygote, ForwardDiff
using HANK.OC

#using ModelingToolkit
using Symbolics

## A script for experimenting and developing
## the method of Mortezaee & Nazemi (2019)
F(z,c) = log(c[1])
r = 0.05
y = 2.0
g = let r = r, y=y
    G(z,c) = [(r*z[1]-c[1] + y)]
end
p = OC.OCProblem(F,g,[0., 0.],[100.,100.], 1,1, [10.0], r)

H = OC.hamiltonian(p)


@variables z[1:p.dimz],c[1:p.dimc]
@variables λ[1:p.dimz],ts
z = Symbolics.scalarize(z)
c = Symbolics.scalarize(c)
λ = Symbolics.scalarize(λ)
∂H∂z = build_function(Symbolics.gradient(H(ts,z,c,λ),z), ts,z,c, λ)[1] |> eval
∂H∂c = build_function(Symbolics.gradient(H(ts,z,c,λ),c), ts,z,c, λ)[1] |> eval  

#h = NonlinearFunction(ns, jac=true)

points = 20
width = 20

import HANK.OC: constraindomain

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

import Random
rng = Random.default_rng()
θc, stc = Lux.setup(rng,cm)
θc.layer_2.bias .= Random.rand(-scalet.(t), width)
θc.layer_3.weight .*= 0.1
θc.layer_3.bias .= -log((p.U[1]-p.L[1])/(1.0+p.r*p.z0[1]-p.L[1])-1) 
θz, stz = Lux.setup(rng,zm)
θz.layer_2.bias .= Random.rand(-scalet.(t), width)
θz.layer_3.bias .= -log((p.U[2]-p.L[2])/(p.z0[1]-p.L[2])-1) 
θλ, stλ = Lux.setup(rng,λm)
θλ.layer_2.bias .= Random.rand(-scalet.(t), width)


focs(t, z, c, λ, ż, λ̇) = vcat(
        ż .- p.G(z,c),
        ∂H∂z(t,z,c,λ) + λ̇,
        ∂H∂c(t,z,c,λ))
t = RType.(t)
import FiniteDiff
function loss(θc, θz, θλ)
    C(t) = cm([t], θc, stc)[1]
    Z(t) = RType.(p.z0) .+ tanh(t)*zm([t], θz, stz)[1]
    Ż(t) = FiniteDiff.finite_difference_derivative(Z,t)
    Λ(t) = exp(-RType(p.r)*t)*λm([t], θλ, stλ)[1]
    #Λ(t) = λm([t], θλ, stλ)[1]
    Λ̇(t) = FiniteDiff.finite_difference_derivative(Λ,t)
    sum( sum(x->x'*x, 
        focs(s, Z(s), C(s), Λ(s), Ż(s), Λ̇(s)))
        for s ∈ t)
end
    

θ = (θc,θz,θλ)
x, re = Lux.destructure(θ)
ℓ = x -> loss(re(x)...)
lval, back = Zygote.pullback(ℓ,x)
@time gs = back(one(lval))[1];
#@show gs
@time gf = ForwardDiff.gradient(ℓ,x);
@show gf ≈ gs

opt = Lux.Optimisers.setup(Lux.Optimisers.ADAM(0.1f0), x)
iterations = 500
for i ∈ 1:iterations
    lval, back = Zygote.pullback(ℓ, x)
    gs = back(one(lval))[1]
    opt, x = Lux.Optimisers.update(opt, x, gs)
    println("Iter[$i]: obj=$(lval)")
end

#x = zero(x)
θc, θz, θλ = re(x)
C(t) = cm([t], θc, stc)[1]
Z(t) = RType.(p.z0) .+ tanh(t)*zm([t], θz, stz)[1]
Ż(t) = FiniteDiff.finite_difference_derivative(Z,t)
Λ(t) = exp(-RType(p.r)*t)*λm([t], θλ, stλ)[1]
Λ̇(t) = FiniteDiff.finite_difference_derivative(Λ,t)
    
using PlotlyLight
plt = Plot()
plt(x=t,y=vcat(C.(t)...), name="consumption", mode="lines+markers")
plt(x=t,y=vcat(Z.(t)...), name="savings", mode="lines+markers")
plt(x=t,y=vcat(Λ.(t)...), name="lambda", mode="lines+markers")
plt

fp = Plot()
err = hcat([vcat(focs(s, Z(s), C(s), Λ(s), Ż(s), Λ̇(s))...) for s ∈ t]...)
fp(x=t, y=err[1,:], name="dz err")
fp(x=t, y=err[2,:], name="c err")
fp(x=t, y=err[3,:], name="dl err")

