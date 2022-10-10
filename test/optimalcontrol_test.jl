@testset "quadratic-consumption-savings" begin
    import HANK: OC
    import InfiniteOpt, Ipopt
    F(z) = log(z[1])
    r = 0.05
    y = 0.0
    g = let r = r, y=y
        G(z,dotz) = dotz[2] - (r*z[2]-z[1] + y)
    end
    p = OC.OCProblem(F,g,[0., 0.],[Inf,Inf], 2, [nothing, 100.0], r)
    sol = OC.solve_InfiniteOpt(p, nsupport = 500)


    using PlotlyLight
    data = Config(t = sol.t, Z = sol.Z, type="scatter", mode="lines+markers")
    plt = Plot(data)
end