import Pkg; Pkg.add("DataStructures")
Pkg.add("MutableArithmetics")
Pkg.add("StarAlgebras")
Pkg.add("PyCall")
include("/Users/ray/Documents/Geo/MultivariatePolynomials.jl/src/MultivariatePolynomials.jl")
using TypedPolynomials
@polyvar a b c x y z
SA = b^2 + c^2 - a^2
SB = a^2 + c^2 - b^2
SC = a^2 + b^2 - c^2
A = [1,0,0]
B = [0,1,0]
C = [0,0,1]
O = [a^2*SA, b^2*SB, c^2*SC]
H = [SB*SC, SA*SC, SA*SB]
I = [a,b,c]
Ia = [-a,b,c]
Ib = [a,-b,c]
Ic = [a,b,c]
Fe = [(b+c-a)*(b-c)^2, (a+c-b)*(a-c)^2, (a+b-c)*(a-b)^2]

function reducepair(v2)
    g = gcd(v2[1], v2[2])
    return [div_multiple(v2[1], g), div_multiple(v2[2], g)]
end

function reduce2(circle)
    cyz = circle[1]
    cxz = circle[2]
    cxy = circle[3]
    cx = circle[4]
    cy = circle[5]
    cz = circle[6]
    g = gcd(gcd(gcd(gcd(gcd(cyz, cy), cx), cxy), cxz), cz)
    return [czy, cxz, cxy, cx, cy, cz]  
    # return [div_multiple(cyz, g),div_multiple(cxz, g),div_multiple(cxy, g),div_multiple(cx, g),div_multiple(cy, g),div_multiple(cz, g)]
end

function reduce(p1)
    polyx = p1[1]
    polyy = p1[2]
    polyz = p1[3]
    g=1
    if(polyx == 0 && polyy == 0)
        return [0,0,1]
    end
    if(polyz == 0 && polyy == 0)
        return [1,0,0]
    end
    if(polyz == 0 && polyx == 0)
        return [0,1,0]
    end
    if(polyx == 0)
        g = gcd(polyy, polyz)
        return [0, div_multiple(polyy,g), div_multiple(polyz,g)]
    end
    if(polyy == 0)
        g = gcd(polyx, polyz)
        return [div_multiple(polyx,g), 0, div_multiple(polyz,g)]
    end
    if(polyz == 0)
        g = gcd(polyy, polyx)
        return [div_multiple(polyx,g), div_multiple(polyy,g), 0]
    end
    g = gcd(polyx, polyy)
    g = gcd(g, polyz)
    return [div_multiple(polyx,g), div_multiple(polyy,g), div_multiple(polyz,g)]
    # return [polyx, polyy, polyz]
end

function midpoint(p1,p2)
    sum1 = p1[3] + p1[1]  + p1[2]
    sum2 = p2[3] + p2[1] + p2[2]
    p3 = [p1[1] * sum2 + p2[1] * sum1, p1[2] * sum2 + p2[2] * sum1, p1[3] * sum2 + p2[3] * sum1]
    p3 = reduce(p3)
    # return p3
end

function intersect(l1,l2)
    x1 = l1[1]
    y1 = l1[2]
    z1 = l1[3]
    x2 = l2[1]
    y2 = l2[2]
    z2 = l2[3]
    pt = [y1*z2 - y2*z1, z1*x2-x1*z2, x1*y2-x2*y1]
    pt = reduce(pt)
    # return pt
end

function area(p1,p2,p3)
    x1 = p1[1]
    y1 = p1[2]
    z1 = p1[3]
    x2 = p2[1]
    y2 = p2[2]
    z2 = p2[3]
    x3 = p3[1]
    y3 = p3[2]
    z3 = p3[3]
    return [x1*y2*z3 + y1*z2*x3 + z1*x2*y3 - x1*z2*y3 - y1*x2*z3 - z1*y2*z3, (x1+y1+z1) * (x2 + y2 + z2) * (x3 + y3 + z3)]
end

function line(p1, p2)
    return intersect(p1,p2)
end

function translate(vec1, p1)
    sum1 = vec1[1]+vec1[2]+vec1[3]
    sum2 = p1[1]+p1[2]+p1[3]
    return reduce([p1[1]*sum1 + vec1[1]*sum2,p1[2]*sum1 + vec1[2]*sum2,p1[3]*sum1 + vec1[3]*sum2])
end

function altitude(l1, p1)
    l2 = perpendicular_line(l1, p1)
    p2 = intersect(l1,l2)
    return reduce(p2)
end

function reflect_point(l1, p1)
    p2 = altitude(l1,p1)
    v = vector(p2, p1)
    return reduce(translate(v, p2))
end

function reflect_line(l1,l2)
    p1 = 
end

function reduce_dist(dist)
    
end

function vector(p1, p2)
    x1 = p1[1]
    y1 = p1[2]
    z1 = p1[3]
    x2 = p2[1]
    y2 = p2[2]
    z2 = p2[3]
    sum1 = x1+y1+z1
    sum2 = x2+y2+z2
    p3 = [x1*sum2 - x2*sum1, y1*sum2 - y2*sum1, z1*sum2 - z2*sum1]
    return p3
end

function perpendicular_line(l1, p1)
    p2 = [-l1[2] - l1[3], l1[1], l1[1]]
    p3 = [l1[2], -l1[1]-l1[3], l1[2]]
    v1 = vector(p2,p3)
    vx = v1[1]
    vy = v1[2]
    vz = v1[3]
    px = p1[1]
    py = p1[2]
    pz = p1[3]
    

    nx = a^2 * (-vy * pz - vz * py) + b^2*(-vx*pz + vz * (py + pz)) + c^2 * (-vx * py + vy *(py+pz))
    ny = a^2 * (-vy * pz + vz * (px+pz)) + b^2*(-vx*pz - vz * px) + c^2 * (vx *(px+pz) - vy *px)
    nz = a^2 * (vy * (py+px) - vz * py) + b^2*(vx*(py+px) - vz * px) + c^2 * (-vx * py - vy *px)
    line = [nx, ny, nz]
    return reduce(line)
    # return line
end

# a^2 b^2 c^2
# vx vy vz 
# s1x - px(x+y+z) s1y - py(x+y+z) s1z - pz(x+y+z)

function circumcenter(p1,p2,p3)
    mid1 = midpoint(p1,p2)    
    mid2 = midpoint(p1,p3)
    p1p2 = line(p1,p2)
    p1p3 = line(p1,p3)
    l1 = perpendicular_line(p1p2, mid1)
    l2 = perpendicular_line(p1p3, mid2)
    center = intersect(l1,l2)
    return center
    # return reduce(center)
end

function sq_distance(p1,p2)
    v = vector(p1,p2) #not reduced
    s = (p1[1] + p1[2] + p1[3]) * (p2[1] + p2[2] + p2[3])
    return [-a^2*v[2]*v[3] - b^2 * v[1] * v[3] - c^2 * v[1] * v[2], s^2]
    # return reducepair([-a^2*v[2]*v[3] - b^2 * v[1] * v[3] - c^2 * v[1] * v[2], s^2])
end

function reduce_circle(circle)
    c1 = circle[1]
    c2 = circle[2]
    c3 = circle[3]
    c4 = circle[4]
    c5 = circle[5]
    c6 = circle[6]
    g = gcd(c1,c2)
    g = gcd(g,c3)
    g = gcd(g,c4)
    g = gcd(g,c5)
    g = gcd(g,c6)
    return [c1,c2,c3,c4,c5,c6]
    # return [div_multiple(c1,g), div_multiple(c2,g),div_multiple(c3,g),div_multiple(c4,g),div_multiple(c5,g),div_multiple(c6,g)]
end

function circle(p1,p2,p3)
    center = circumcenter(p1,p2,p3)
    dist = sq_distance(center, p1)
    cx = center[1]
    cy = center[2]
    cz = center[3]
    sum = cx+cy+cz
    ryz = -a^2 * dist[2] * sum
    rxz = -b^2 * dist[2] * sum
    rxy = -c^2 * dist[2] * sum
    rx = b^2 * cz * dist[2] * sum + dist[2] * c^2 * cy * sum - dist[1] * sum^2 - (a^2 * cz * cy + b^2 * cx * cz + c^2 * cx * cy)*dist[2]
    ry = a^2 * cz * dist[2] * sum + dist[2] * c^2 * cx * sum - dist[1] * sum^2 - (a^2 * cz * cy + b^2 * cx * cz + c^2 * cx * cy)*dist[2]
    rz = a^2 * cy * dist[2] * sum + dist[2] * b^2 * cx * sum - dist[1] * sum^2 - (a^2 * cz * cy + b^2 * cx * cz + c^2 * cx * cy)*dist[2]
    # return reduce_circle([ryz, rxz, rxy, rx, ry, rz, sum])
    return [ryz, rxz, rxy, rx, ry, rz, sum]
end

# sum * (-a^2(y-j)(z-k) - b^2(x-i)(z-k) - c^2(x-i)(y-j)) = dist * sum^2 (x+y+z) # x+y+z = sum

function intersect_lc(circ, l)
    # sum = circ[7]
    # x + y + z = sum, x = sum- y - z
    l1 = x*l[1] + y * l[2] + z*l[3]
    print(l1)
end

# x+y+z = sum
# px+qy+rz = 0
# intersect
# c1 yz + c2 xz + c3 xy + c4 x + c5 y + c6 z = 0

function intersect_cc(circ1, circ2)

end