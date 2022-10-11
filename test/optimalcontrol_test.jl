@testset "quadratic-consumption-savings" begin
    import HANK: OC
    import InfiniteOpt, Ipopt, Lux
    F(z) = log(z[1])
    r = 0.05
    y = 0.0
    g = let r = r, y=y, α=0.6
        G(z,dotz) = dotz[2] - (r*z[2]-z[1] + y)
    end
    p = OC.OCProblem(F,g,[0., 0.],[Inf,Inf], 2, [1.0, 1.0], r)
    sol = OC.solve_InfiniteOpt(p, nsupport = 1000)


    using PlotlyLight
    plt = Plot()
    plt(x=sol.t,y=sol.Z[1,:], name="consumption", mode="lines+markers")
    plt(x=sol.t,y=sol.Z[2,:], name="savings", mode="lines+markers")
    plt

    m = OC.controlmodel(p; width=10)
    θ, st = Lux.setup(Random.default_rng(), m)
    θa, re = Lux.destructure(θ)
    ode = OC.createode(p, m, st)(θa) 
    s = OC.solveode(ode,θa)
    plt = Plot()
    plt(x=sol.t,y=vcat((t->m([t], re(θa), st)[1][1]).(sol.t)...), name="consumption", mode="lines+markers")
    plt(x=sol.t,y=vcat(s.(sol.t)...), name="savings", mode="lines+markers")
    plt

    OC.traincontrol!(θ, st, m, p, nint=10, iterations=10)
       
    plt = Plot()
    plt(x=sol.t,y=vcat((t->m([t], θ, st)[1][1]).(sol.t)...), name="consumption", mode="lines+markers")
    θa, re = Lux.destructure(θ)
    s = OC.solveode(ode,θa)
    plt(x=sol.t,y=vcat(s.(sol.t)...), name="savings", mode="lines+markers")
    plt
    

end