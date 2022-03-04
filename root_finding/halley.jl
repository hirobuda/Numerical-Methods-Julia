function halley(f::Function, f1::Function, f2::Function, x0::Float64)
    
    MAX = 1000
    x = zeros(MAX)
    x[1] = x0
    
    for i=2:MAX
        x[i] = x[i-1] - (f(x[i-1])*f1(x[i-1]))/(f1(x[i-1])^2 - 0.5*f(x[i-1])*f2(x[i-1])) 
        if abs(x[i] - x[i-1]) < 10^(-16)
            return x
        end        
    end
    
    return x
end
