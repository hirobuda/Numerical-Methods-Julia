#u_t + u_x = 0
#u(x,0) = η(x) = exp(-20(x-2)^2) + exp(-(x-5)^2)
function adv_frog(h,k,xf,tf)
    
    m = Int(ceil(xf/h))
    n = Int(ceil(tf/k))
    r = k/h
    
    A = Tridiagonal(r*ones(m-1), zeros(m), -r*ones(m-1))
    
    x = LinRange(0,xf,m)
    t = zeros(n)
    t[1] = k
    
    U = zeros(m,n)
    U[:,1] = η.(x) 
    A_lax = Tridiagonal((0.5+r)*ones(m-1), zeros(m), (0.5-r)*ones(m-1))
    U[:,2] = A_lax*U[:,1]
    
    for i=3:n
        U[:,i] = A*U[:,i-1] + U[:,i-2]
        t[i] = i*k
    end
    
    return m,n,x,t,U
end
