function secant(f::Function, x0::Float64, x1::Float64)::Array{Float64}
    
    MAX = 1000
    x = zeros(MAX)
    x[1] = x0
    x[2] = x1
    
    for i=3:MAX
        x[i] = (x[i-2]*f(x[i-1]) - x[i-1]*f(x[i-2]))/(f(x[i-1]) - f(x[i-2]))
        if abs(x[i] - x[i-1]) < 10^(-16)
            return x
        end
    end
    
    return x
end
