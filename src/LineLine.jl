function LineLine(l1::Geometry, l2::Geometry)
    #Ips = Array{Tuple{Point2D{Float64}, Float64, Float64}, 1}()
    dx_cx = l2.Ep.x - l2.Sp.x
    dy_cy = l2.Ep.y - l2.Sp.y
    bx_ax = l1.Ep.x - l1.Sp.x
    by_ay = l1.Ep.y - l1.Sp.y
    Δ = bx_ax * dy_cy - by_ay * dx_cx
    if abs(Δ) < 1e-4 # colinear
        if l1.Ep.x == l2.Sp.x # l2 follows l1
            return (intersect=true, type="connect", Ip=Point4To2D(l1.Ep), Dist1=l1.Length, Dist2=0.0)
        elseif l2.Ep.x == l1.Sp.x # l1 follows l2
            return (intersect=true, type="connect", Ip=Point4To2D(l1.Sp), Dist1=0.0, Dist2=l2.Length)
        elseif onSegment(l1.Sp, l1.Ep, l2.Sp)
            return (intersect=true, type="overlap", Ip=Point4To2D(l2.Sp), Dist1=Distance2D(l1.Sp, l2.Sp), Dist2=0.0)
        elseif onSegment(l2.Sp, l2.Ep, l1.Sp)
            return (intersect=true, type="overlap", Ip=Point4To2D(l1.Sp), Dist1=0.0, Dist2=Distance2D(l2.Sp, l1.Sp))
        else
            return (intersect=false, type="none")
        end
    else # not colinear
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
            return (intersect=true, type="cross", Ip=Ipc, Dist1=dist1, Dist2=dist2)
        else
            return (intersect=false, type="none")
        end
    end
end

function onSegment(Sp::Point4D{Float64}, Ep::Point4D{Float64}, Pi::Point4D{Float64})
    if (Pi.x <= max(Sp.x, Ep.x) && Pi.x >= min(Sp.x, Ep.x) &&
        Pi.y <= max(Sp.y, Ep.y) && Pi.y >= min(Sp.y, Ep.y))
        return true
    else
        return false
    end
end
