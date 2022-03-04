#u_t + u_x = 0
#u(x,0) = η(x) = exp(-20(x-2)^2) + exp(-(x-5)^2)
function adv_lax_wendroff(h,k,xf,tf)
    
    m = Int(ceil(xf/h))
    n = Int(ceil(tf/k))
    r = k/h
    
    A = Tridiagonal(0.5*(r^2+r)*ones(m-1), (1-r^2)*ones(m), 0.5*(r^2-r)*ones(m-1))
    
    x = LinRange(0,xf,m)
    t = zeros(n)
    t[1] = k
    t[2] = 2*k
    U = zeros(m,n)
    U[:,1] = η.(x) 
    for i=2:n
        U[:,i] = A*U[:,i-1] 
        t[i] = i*k
    end
    
    return m,n,x,t,U
end
