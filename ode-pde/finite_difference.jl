# u''(x) = exp(x), 0 <= x <= 1
# u'(0) = 0, u(1) = 3

function finite_difference(h::Float64, sigma::Float64, beta::Float64, f::Function)
    
    m = Int(ceil(1/h))

    d_up = ones(m+1)
    d_up[1] = h

    d = -2*ones(m+2)
    d[1] = -h
    d[m+2] = h^2

    d_down = ones(m+1)
    d_down[m+1] = 0

    A = Tridiagonal(d_down,d,d_up)
    A = (1/h^2)*A
    
    x = LinRange(0,1,m+2)
    
    G = zeros(m+2)
    G[1:m+2] = f.(x)
    G[1] = sigma + (h/2)*f(x[1])
    G[m+2] = beta
    
    return A\G, x
end
