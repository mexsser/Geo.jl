struct Point2D{T}
    x::T
    y::T
end

Base.:-(P1::Point2D{T}, P2::Point2D{T}) where T = Point2D{T}(P1.x-P2.x, P1.y-P2.y)
Base.:+(P1::Point2D{T}, P2::Point2D{T}) where T = Point2D{T}(P1.x-P2.x, P1.y-P2.y)

struct Point4D{T}
    x::T
    y::T
    cur::T  # Curvature at this point, positive if dθ/ds > 0
    θ::T  # tangential direction in global coordinates
end

function Point4To2D(P::Point4D{T}) where T
    return Point2D{T}(P.x, P.y)
end

function Distance2D(Sp::Point2D{Float64}, Ep::Point2D{Float64})
    return sqrt((Sp.x-Ep.x)^2 + (Sp.y-Ep.y)^2)
end

function Distance2D(Sp::Array{Float64, 1}, Ep::Array{Float64, 1})
    return sqrt((Sp[1]-Ep[1])^2 + (Sp[2]-Ep[2])^2)
end

function Distance2D(Sp::Point4D{Float64}, Ep::Point4D{Float64})
    return sqrt((Sp.x-Ep.x)^2 + (Sp.y-Ep.y)^2)
end

function round(x::Float64, pre::Union{Float64, Int64}) # pre: precision
    n = floor(Int, x/pre)
    m = x - n*pre
    if m >= pre/2
        return trunc((n+1)*pre; digits=10)
    else
        return trunc(n*pre; digits=10)
    end
end

struct Geometry
    Sp::Point4D{Float64} # Start Point4D
    Ep::Point4D{Float64} # End Point4D
    Length::Float64
    Vstd::Float64 # the standard velocity of this Geometry
    GType::String
end

function Geometry(Sp::Point4D{Float64}, Ep::Point4D{Float64}, Vstd::Float64)
    if Sp.cur == Ep.cur
        if Sp.cur != 0.0
            Length = round(abs((Sp.θ - Ep.θ)/Sp.cur), 1.0e-14)
            GType = "arc"
        else
            return Geometry(Point2D{Float64}(Sp.x, Sp.y), Point2D{Float64}(Ep.x, Ep.y), Vstd) #GType = "line"
        end
    else
        GType = "clothoid"
        Length = 1.0
    end
    return Geometry(Sp, Ep, Length, Vstd, GType)
end

function Geometry(Sp::Point2D{Float64}, Ep::Point2D{Float64}, Vstd::Float64) # constructor for line element
    Length = round(Distance2D(Sp, Ep), 1.0e-14)
    Δx = Ep.x-Sp.x
    Δy = Ep.y-Sp.y
    Theta = 0.0
    if Δx != 0
        Theta = atan(Δy/Δx)
    elseif Δy > 0
        Theta = π/2
    elseif Δy < 0
        Theta = -π/2
    else
        error("Start point and end point are the same")
    end
    return Geometry(Point4D{Float64}(Sp.x, Sp.y, 0.0, Theta), Point4D{Float64}(Ep.x, Ep.y, 0.0, Theta), Length, Vstd, "line") # using default constructor
end

function ArcCenter(arc::Geometry)
    if arc.GType != "arc"
        error("input type incorrect: arc needed.")
    end
    radius = abs(1/arc.Sp.cur)
    if arc.Sp.θ > arc.Ep.θ
        δ = π/2
    else
        δ = -π/2
    end
    Sp2 = [radius; 0]
    Rx = [cos(arc.Sp.θ + δ) -sin(arc.Sp.θ + δ); sin(arc.Sp.θ + δ) cos(arc.Sp.θ + δ)]
    Result = [arc.Sp.x; arc.Sp.y] - Rx*Sp2
    Cpoint = Point2D{Float64}(Result[1], Result[2])
    return Cpoint
end

function PointInSector(Ip::Point2D{Float64}, arc::Geometry)
    Cpoint = ArcCenter(arc)
    Vip = [Ip.x-Cpoint.x, Ip.y-Cpoint.y]   # due to matrix malipulation, self-defiend types should not be used
    Vsp = [arc.Sp.x-Cpoint.x, arc.Sp.y-Cpoint.y]
    Vep = [arc.Ep.x-Cpoint.x, arc.Ep.y-Cpoint.y]
    if Distance2D(Vip, Vsp) <= 1.0e-6
        return (true, 0.0)
    elseif Distance2D(Vip, Vep) <= 1.0e-6
        return (true, arc.Length)
    else
        φs = atan(Vsp[2], Vsp[1]) # [-π, π]
        Rx = [cos(φs) sin(φs); -sin(φs) cos(φs)]
        Vip_r = Rx*Vip  # rotate
        Vep_r = Rx*Vep
        φe_r = atan(Vep_r[2], Vep_r[1])
        φi_r = atan(Vip_r[2], Vip_r[1])

        if arc.Sp.cur > 0 # from Sp to Ep: counter-clockwise
            if (φe_r > 0 && 0 < φi_r < φe_r) || (φe_r < 0 && 0 < φi_r < φe_r+2π)
                return (true, φi_r/arc.Sp.cur)
            end
        else # from Sp to Ep: clockwise
            if (φe_r < 0 && φe_r < φi_r < 0) || (φe_r > 0 && φe_r-2π < φi_r < 0)
                return (true, φi_r/arc.Sp.cur)
            end
        end
    end
    return (false, nothing)
end

function Length2XY(Geo::Geometry, Len::Float64)  # convert local 1D coordinate(:length) to global 2D coordinate. The heading direction θ can also be returned if needed
    Δx = 0.0
    Δy = 0.0

    if Geo.GType == "arc"
        radius = 1/Geo.Sp.cur # no need to calculate radius using (x,y,θ), just get it from curvature
        ΔΘ = Len/radius
        Center = ArcCenter(Geo)
        RotMat = [cos(ΔΘ) -sin(ΔΘ); sin(ΔΘ) cos(ΔΘ)]
        Δx, Δy = ([-1 0; 0 -1]+RotMat)*[Geo.Sp.x - Center.x; Geo.Sp.y - Center.y]
        #Δx = sign*radius*( sin(Geo.Sp.θ+sign*ΔΘi) - sin(Geo.Sp.θ) )
        #Δy = sign*radius*( cos(Geo.Sp.θ) - cos(Geo.Sp.θ+sign*ΔΘi) )
    elseif Geo.GType == "line" # straight line
        Δx = (Geo.Ep.x - Geo.Sp.x) * Len/Geo.Length
        Δy = (Geo.Ep.y - Geo.Sp.y) * Len/Geo.Length
    else # clothoid
        # construct later
    end
    return Point2D{Float64}(Geo.Sp.x+Δx, Geo.Sp.y+Δy)
end

function Length2XYθ(Geo::Geometry, Len::Float64)
    Point = Length2XY(Geo, Len)
    θ = Geo.GType == "arc" ? Geo.Sp.θ + Len*Geo.Sp.cur : Geo.Sp.θ
    return Point, θ
end

# the following code block should be improved using better strcuts' structure
function GeosIntersect(Geo1::Geometry, Geo2::Geometry)    # consider using enumerate to handle names
    if Geo1.GType == Geo2.GType
        if Geo1.GType == "line"
            return LineLine(Geo1, Geo2)
        elseif Geo1.GType == "arc"
            return ArcArc(Geo1, Geo2)
        end
    else
        if Geo1.GType == "line" && Geo2.GType == "arc"
            return LineArc(Geo1, Geo2)
        elseif Geo2.GType == "line" && Geo1.GType == "arc"
            return LineArc(Geo2, Geo1)
        end
    end
end
