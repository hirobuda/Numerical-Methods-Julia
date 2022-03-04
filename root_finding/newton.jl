function newton(f::Function, df::Function, x0::Float64)
    
    MAX = 1000
    x = zeros(MAX)
    x[1] = x0
    
    for i=2:MAX
        x[i] = x[i-1] -f(x[i-1])/df(x[i-1])
        if abs(x[i] - x[i-1]) < 10^(-16)
            return x
        end        
    end
    
    return x
end
