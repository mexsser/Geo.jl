function LineLine(l1::Geometry, l2::Geometry)
    #Ips = Array{Tuple{Point2D{Float64}, Float64, Float64}, 1}()
    dx_cx = l2.Ep.x - l2.Sp.x
    dy_cy = l2.Ep.y - l2.Sp.y
    bx_ax = l1.Ep.x - l1.Sp.x
    by_ay = l1.Ep.y - l1.Sp.y
    Δ = bx_ax * dy_cy - by_ay * dx_cx
    if abs(Δ) < 1e-4 # parallel
        if l1.Ep.x == l2.Sp.x # l2 follows l1
            return (intersect=true, Ip=Point4To2D(l1.Ep), Dist1=l1.Length, Dist2=0.0)
        elseif l2.Ep.x == l1.Sp.x # l1 follows l2
            return (intersect=true, Ip=Point4To2D(l1.Sp), Dist1=0.0, Dist2=l2.Length)
        else
            return (intersect=false, message="line overlap")
        end
    else # not paralell
        ax_cx = l1.Sp.x - l2.Sp.x
        ay_cy = l1.Sp.y - l2.Sp.y
        r = (ay_cy * dx_cx - ax_cx * dy_cy) / Δ
        s = (ay_cy * bx_ax - ax_cx * by_ay) / Δ
        if (0.0 <= r <= 1.0 &&  0.0 <= s <= 1.0)  # intersection within 2 line-segments
            Ix = l1.Sp.x + r * bx_ax
            Iy = l1.Sp.y + r * by_ay
            Ipc = Point2D(Ix, Iy)
            dist1 = Distance2D(Point4To2D(l1.Sp), Ipc)
            dist2 = Distance2D(Point4To2D(l2.Sp), Ipc)
            return (intersect=true, Ip=Ipc, Dist1=dist1, Dist2=dist2)
        else
            return (intersect=false, message="nothing")
        end
    end
    #push!(Ips, Ip)
    #return (true, Ips)
end
