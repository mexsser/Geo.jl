function LineArc(line::Geometry, arc::Geometry)
    if line.GType != "line" || arc.GType != "arc"
        error("input type incorrect: line and arc needed")
    end
    Center = ArcCenter(arc)
    SymPy.@vars x y t
    x = line.Sp.x + t*(line.Ep.x - line.Sp.x)
    y = line.Sp.y + t*(line.Ep.y - line.Sp.y)
    Eq = (x - Center.x)^2 + (y - Center.y)^2 - (1/arc.Sp.cur)^2
    Roots = SymPy.solve(Eq)
    #Ips = Array{Tuple{Point2D{Float64},Float64, Float64}, 1}()
    if !isa(SymPy.N(Roots[1]), Real)  # line and circle do not Intersect
        #println(typeof(Roots[1]))
        return (intersect=false, message="no real roots")
    else
        for ti in Roots
            if 0.0 <= ti <= 1.0 + 1.0e-6  # check if intersection point is in line segment
                Ipc = Point2D{Float64}(SymPy.subs(x, t, ti), SymPy.subs(y, t, ti))
                dist1 = SymPy.N(ti*line.Length)
                Result = PointInSector(Ipc, arc)
                if Result[1] # check if intersection point is in arc
                    dist2 = Result[2]
                    return (intersect=true, Ip=Ipc, Dist1=dist1, Dist2=dist2)
                    #push!(Ips, (Ipc, Dist1, Dist2))
                end
            end
        end
    end

    return (intersect=false, message="no intersection")
#    if size(Ips, 1) == 0
#        return (false, "not inside")
#    else
#        return (true, Ips)
#    end
end
