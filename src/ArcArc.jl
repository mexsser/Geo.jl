function ArcArc(arc1::Geometry, arc2::Geometry)
    if arc1.GType != "arc" || arc2.GType != "arc"
        error("Input Type error： arc needed")
    end
    Cp1 = ArcCenter(arc1)
    Cp2 = ArcCenter(arc2)
    r1 = 1/abs(arc1.Sp.cur)
    r2 = 1/abs(arc2.Sp.cur)
    Dis = Distance2D(Cp1, Cp2)
    if Dis > r1 + r2 || Dis < abs(r1 - r2) # separate or contain
        return (intersect=false, type="none")
    elseif Dis <= 1.0e-6 && r1 == r2 # coincident || overlap
        error("arcs are in the same circle, may be handled later.")
    else
        a = (r1^2 - r2^2 + Dis^2)/(2*Dis)
        h = sqrt(r1^2 - a^2)
        xm = Cp1.x + a*(Cp2.x - Cp1.x)/Dis
        ym = Cp1.y + a*(Cp2.y - Cp1.y)/Dis
        Ipx1 = xm + h*(Cp2.y - Cp1.y)/Dis
        Ipy1 = ym - h*(Cp2.x - Cp1.x)/Dis
        Ipx2 = xm - h*(Cp2.y - Cp1.y)/Dis
        Ipy2 = ym + h*(Cp2.x - Cp1.x)/Dis
        Ip1 = Point2D{Float64}(Ipx1, Ipy1)
        Ip2 = Point2D{Float64}(Ipx2, Ipy2)
        Infos = Array{NamedTuple{(:Ip, :Dist1, :Dist2),Tuple{Point2D{Float64}, Float64, Float64}}, 1}()
        for Ipi in [Ip1, Ip2]
            Result1 = PointInSector(Ipi, arc1)
            Result2 = PointInSector(Ipi, arc2)
            if Result1[1] && Result2[1]
                push!(Infos, (Ip=Ipi, Dist1=Result1[2], Dist2=Result2[2]))
            end
        end
        if length(Infos) == 0
            return (intersect=false, type="none")
        else
            return (intersect=true, type="cross", Infos[1]...) # actually both points should be used to calculate TTC. But now I just use one and this will be changed later
        end
    end
end
