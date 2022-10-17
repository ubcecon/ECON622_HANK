##############################################
############# using packages #################
##############################################
using NLsolve
using Plots
using LinearAlgebra
using EconPDEs
using JLD
import Random
################################################
############# set parameters ###################
################################################

w = 1.0 
T = 0.0 
r_a = 0.07
r_b = 0.05 
β = 0.02
ρ = 0.051
λ = 0.0 
γ = 2.0 
v = 1.0
χ_0 = 0.0438
χ_1 = 0.956
χ_2 = 1.402
ϕ = 2.2

function initialize_y(stategrid)
    an = length(stategrid[:a])
    bn = length(stategrid[:b])
    yn = length(stategrid[:y])
    OrderedDict(:V => ones(an, bn, yn))
end



as = collect(range(0.0, stop=10.0, length=100))
bs = collect(range(-10.0, stop=10.0, length=200))
ys = collect(range(0, stop=2, length=100))
stategrid = OrderedDict(:a => as, :b => bs, :y => ys)
yend = initialize_y(stategrid)

function f(state::NamedTuple, sol::NamedTuple)
    a= state.a 
    b = state.b
    y = state.y
    V, Va_up, Va_down, Vb_up, Vb_down, Vy_up, Vy_down = sol.V, sol.Va_up, sol.Va_down, sol.Vb_up, sol.Vb_down, sol.Vy_up, sol.Vy_down 

    # upwinding
    Va = Va_up
    Vb = Vb_up
    Vy = Vy_up

    # conpute
    c = (Vb)^(-1/γ)
    l = ((1/ϕ)* Vb* (1-τ) *w * (e^y))^(1/v)

    d = ((Va /Vb - 1) - χ_0)/(χ_1*χ_2* a^(1-χ_2))^(1/(χ_2-1)) > 0 ? ((Va /Vb - 1) - χ_0)/(χ_1*χ_2* a^(1-χ_2))^(1/(χ_2-1)) : ((Va /Vb - 1) * (-1) - χ_0)/(χ_1*χ_2* a^(1-χ_2))^(1/(χ_2-1))


    # define a function 
    χ = χ_0 * abs(d) + χ_1*((abs(d/a))^χ_2)*a
    # solving endogeneous parameters 
    Vt = - ( (c^(1-γ))/(1-γ) - ϕ * ((l^(1+v))/(1+v)) + Vb * ((1-τ)*w* exp(y) *l + r_b * b +T - d - χ -c) + Va * (r_a*a + d) + Vy*(-β*y) - (ρ + ζ)*V )
    (Vt = Vt,), (; V, Va, Vb, Vy)
end




result = pdesolve(f, stategrid, yend)

(V = result.optional[:V], Va =result.optional[:Va], Vb = result.optional[:Vb], Vy = result.optional[:Vy])    

# Phillips Curve
# u = [r_a(t) Y(t) P(t) π(t) m(t) r_b(t)] vector of functions
# z = [ϵ θ rbfixed ϕ] vector of parameters
function PhillipsCurve(du,u,z,t)
    (u[1]- (du[2]/u[2]))*u[4] = (z[1]/z[2])*(u[5]- (z[1]-1)/z[1]) + du[4]
    du[3]/u[3] = u[4]
    u[6] = z[3] + (z[4]-1)*u[4] #Taylor rule

   end
#Goverment's Law of Motion 
# u = [Bg(t) G(t) T(t) τ(t) w(t) l(t) rb(t)] vector of functions

function Goverment(du,u,z,t)
 du[1] + u[2] + u[3] = u[4]*∫(u[5]*z*u[6]) + u[7]*u[1]
 end 