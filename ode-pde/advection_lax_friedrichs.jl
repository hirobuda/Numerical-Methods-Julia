function adv_lax_friedrichs(h,k,xf,tf)
    
    m = Int(ceil(xf/h))
    n = Int(ceil(tf/k))
    r = k/h
    
    A = Tridiagonal((0.5+0.5*r)*ones(m-1), zeros(m), (0.5-0.5*r)*ones(m-1))
    
    x = LinRange(0,xf,m)
    t = zeros(n)
    t[1] = k
    
    U = zeros(m,n)
    U[:,1] = η.(x) 
    for i=2:n
        U[:,i] = A*U[:,i-1] 
        t[i] = i*k
    end
    
    return m,n,x,t,U
end
