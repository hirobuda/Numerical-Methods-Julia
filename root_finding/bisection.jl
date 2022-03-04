function bisection(f::Function, a0::Float64, b0::Float64)
    
    MAX = 1000
    x = zeros(MAX)
    a = a0
    b = b0
    
    for i=1:MAX
        c = (a + b)/2
        if abs(a-b) < 10^(-16)
            break        
        elseif f(a)*f(c) < 0
            x[i] = c
            b = c
        elseif f(b)*f(c) < 0
            x[i] = c
            a = c
        end
    end
    
    return x
end

