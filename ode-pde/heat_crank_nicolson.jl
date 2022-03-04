function heat_cn(h,k,n)
    
    r = k/(2*h^2)
    m = Int(ceil(1/h))+1
    x = LinRange(0,1,m)
    A = Tridiagonal(-r*ones(m-3), (1+2*r)*ones(m-2), -r*ones(m-3))
    B = Tridiagonal(r*ones(m-3), (1-2*r)*ones(m-2), r*ones(m-3))
    function f(x)
        if x<=0.5
            return 2*x
        end
        
        if x > 0.5
            return 2-2*x
        end
    end
    U = zeros(m,n)
    U[:,1] .= f.(x)
    t = zeros(n)
    t[1] = 0 
    for j=2:n
        b = zeros(m-2)
        b[1] = r*U[1,j-1]
        b[m-2] = r*U[m,j-1]
        rhs = B*U[2:m-1,j-1] + b
        
        U[2:m-1,j] = A\rhs
        
        t[j] = j*k
    end
    
    return m,n,x,t,U
end
