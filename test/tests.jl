using Geo

function main()
    L1 = Geometry(Point2D{Float64}(1, 0), Point2D{Float64}(1, 2), 4.0)
    L2 = Geometry(Point2D{Float64}(-1, 0), Point2D{Float64}(2, 0), 4.0)
    println("#### Line1 and Line2:")
    println(GeosIntersect(L1, L2))

    L3 = Geometry(Point2D{Float64}(0.0, -1.0), Point2D{Float64}(0.5, 0.0), 4.0)
    A1 = Geometry(Point4D{Float64}(0.0, 0.0, 0.5, -π/6), Point4D{Float64}(2.0, 0.0, 0.5, π/6), 3.0)
    println("\n#### Line3 and Arc1:")
    println(GeosIntersect(L3, A1))

    A2 = Geometry(Point4D{Float64}(0.0, 0.0, 0.5, -π/6), Point4D{Float64}(2.0, 0.0, 0.5, π/6), 3.0)
    A3 = Geometry(Point4D{Float64}(-2.0, 0.0, 0.5, -π/6), Point4D{Float64}(0.0, 0.0, 0.5, π/6), 3.0)
    println("\n#### Arc2 and Arc3:")
    println(GeosIntersect(A2, A3))
end

main()
