include("IRK8.jl")

function ThreeBodyEnergy(u)
    norm12 = sqrt(u[1]^2+u[2]^2)
    norm23 = sqrt(u[3]^2+u[4]^2)
    norm31 = sqrt(u[5]^2+u[6]^2)
    return (u[7]^2+u[8]^2+u[9]^2+u[10]^2+u[11]^2+u[12]^2)/2 - 
           (1/norm12 + 1/norm23 + 1/norm31)
end

function ThreeBodyODE!(du,u,p,t)
    q12 = [u[1], u[2]]
    q23 = [u[3], u[4]]
    q31 = [u[5], u[6]]
    aux12 = q12*dot(q12, q12)^-1.5
    aux23 = q23*dot(q23, q23)^-1.5
    aux31 = q31*dot(q31, q31)^-1.5
    dv1 = aux12 - aux31
    dv2 = aux23 - aux12 
    dv3 = aux31 - aux23       
    @. du[1:4] = u[9:12] - u[7:10] 
    @. du[5:6] = u[7:8] - u[11:12]
    du[7:12] .= [dv1[1], dv1[2], dv2[1], dv2[2], dv3[1], dv3[2]]
    du[13] = 1
    return nothing
end

function ThreeBodyODEGlobalTR!(du,u,p,t)
    q12 = [u[1], u[2]]
    q23 = [u[3], u[4]]
    q31 = [u[5], u[6]]
    dot12 = dot(q12, q12)
    dot23 = dot(q23, q23)
    dot31 = dot(q31, q31)
    norm12 = sqrt(dot12)
    norm23 = sqrt(dot23)
    norm31 = sqrt(dot31)
    aux12 = q12/(dot12*norm12)
    aux23 = q23/(dot23*norm23)
    aux31 = q31/(dot31*norm31)
    dv1 = aux12 - aux31
    dv2 = aux23 - aux12
    dv3 = aux31 - aux23    
    @. du[1:4] = u[9:12] - u[7:10] 
    @. du[5:6] = u[7:8] - u[11:12]
    du[7:12] .= [dv1[1], dv1[2], dv2[1], dv2[2], dv3[1], dv3[2]]
    A = 1/dot12+1/dot23+1/dot31
    B = 1/norm12+1/norm23+1/norm31
    s = (A*B)^-0.5
    @. du *= s  
    du[13] = s  
    return nothing
end

abs2rel(u) = [u[1]-u[3], u[2]-u[4], u[3]-u[5], u[4]-u[6], u[5]-u[1], u[6]-u[2],
                 u[7]-u[9], u[8]-u[10], u[9]-u[11], u[10]-u[12], u[11]-u[7], u[12]-u[8]]

function rel2abs(u)
    x1 = (u[5]-u[1])/3.
    x2 = (u[1]-u[3])/3.
    x3 = (u[3]-u[5])/3.
    y1 = (u[6]-u[2])/3.
    y2 = (u[2]-u[4])/3.
    y3 = (u[4]-u[6])/3.
    
    vx1 = (u[11]-u[7])/3.
    vx2 = (u[7]-u[9])/3.
    vx3 = (u[9]-u[11])/3.
    vy1 = (u[12]-u[8])/3.
    vy2 = (u[8]-u[10])/3.
    vy3 = (u[10]-u[12])/3.
    
    return [x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3]
end

function visualize(u0, odef, p, T0, Tend, n, m, title="", triangles=false)
    odef = ThreeBodyODEGlobalTR!
    p = nothing
    n = 256
    m = 16;
    sol = IRK8(u0, T0, Tend, n, m, odef, p);
    uu = sol.u;
    pl = plot(title=title, aspect_ratio=1)

    uu_abs = [rel2abs(u) for u in uu]
    for j in 1:3
        xx = [u[2*(j-1)+1] for u in uu_abs]
        yy = [u[2*(j-1)+2] for u in uu_abs]

        plot!(xx,yy,legend=false)
    end
    
    if triangles
        u0_abs = uu_abs[1]
        uT_abs = uu_abs[end]
        plot!([u0_abs[i] for i in [1, 3, 5, 1]], [u0_abs[i] for i in [2, 4, 6, 2]], color="gray")
        plot!([uT_abs[i] for i in [1, 3, 5, 1]], [uT_abs[i] for i in [2, 4, 6, 2]], color="purple")
    end

    display(pl)
    #return sol.retcode
end