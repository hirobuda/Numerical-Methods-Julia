# u_t = u_xx, 0 < x < 1
# u(0,t) = u(1,t) = 0, t > 0
# u(x,0) = 2x, if 0 <= x <= 0.5
#          2 - 2x, if 0.5 <= x <= 1
# exact solution is (8/pi^2)*sum([(1/i^2)*sin(0.5*i*pi)*sin(i*pi*x)*exp(-(i^2)*(pi^2)*t) for i=1:infinity])

function heat_euler(h,k,n)
    r = k/(h^2)
    m = Int(ceil(1/h))+1
    x = LinRange(0,1,m)
    A = Tridiagonal(r*ones(m-3), (1-2*r)*ones(m-2), r*ones(m-3))
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
    t[1] = k
    for j=2:n
        b = zeros(m-2)
        b[1] = r*U[1,j-1]
        b[m-2] = r*U[m,j-1]
        U[2:m-1,j] = A*U[2:m-1,j-1] + b
        
        t[j] = j*k
    end
    
    return m,n,x,t,U
end
