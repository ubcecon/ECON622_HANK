@testset "Example tests" begin
    for x ∈ [1, 0.0, rand(), 1//4]
        @test totest(x) ≈ x
    end
    for x ∈ [rand(10), rand(2,3), (1,3)]
        @test totest(x) ≈ sum(x)
    end
end

@testset "Broken" begin
    @test totest([1]) ≈ 2
end