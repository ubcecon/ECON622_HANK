@testset "quadratic-consumption-savings" begin
    import HANK: OC
    import InfiniteOpt, Ipopt, Lux
    F(z,c) = log(c[1])
    r = 0.05
    y = 1.0
    g = let r = r, y=y
        G(z,c) = [(r*z[1]-c[1] + y)]
    end
    p = OC.OCProblem(F,g,[0., 0.],[100.,100.], 1,1, [10.0], r)
    sol = OC.solve_InfiniteOpt(p, nsupport = 1000)


    using PlotlyLight
    plt = Plot()
    plt(x=sol.t,y=sol.C[1,:], name="consumption", mode="lines+markers")
    plt(x=sol.t,y=sol.Z[1,:], name="savings", mode="lines+markers")
    plt

    import Random, Lux
    m = OC.controlmodel(p; width=5)
    θ, st = Lux.setup(Random.default_rng(), m)
    ode = OC.createode(p,θ, m, st)
    θa, re = Lux.destructure(θ)
    θz = zero(θa)
    θz[1:(end-1)] = randn(length(θz)-1)*0.1
    θz[end] = -log((p.U[1]-p.L[1])/(y+p.r*p.z0[1]-p.L[1])-1)
    s = OC.solveode(ode(θz))
    plt = Plot()
    t = range(0.01, 90.0, 1000)
    plt(x=t,y=vcat((t->m([t], re(θz), st)[1][1]).(t)...), name="consumption", mode="lines+markers")
    plt(x=t,y=vcat(s.(t)...), name="savings", mode="lines+markers")    
    plt

    θ = re(θz)
    
    θ, m, st, opt = OC.traincontrol(θ, st, m, p, nint=20, iterations=20, δ = 0.1, μ=0.01)
       
    plt = Plot()
    plt(x=t,y=vcat((t->m([t], θ, st)[1][1]).(t)...), name="consumption", mode="lines+markers")
    θa, re = Lux.destructure(θ)
    s = OC.solveode(ode(θa))
    plt(x=t,y=vcat(s.(t)...), name="savings", mode="lines+markers")
    plt
    

end