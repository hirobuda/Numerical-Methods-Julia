# Lotka-Volterra Equations

F(U,t) = [-2*U[1]+U[1]*U[2], U[2]-U[1]*U[2]]
function predictor_corrector(F,t0,x0,y0,h,N)
    t = zeros(N)
    U = zeros(N,2) 
    
    t[1] = t0
    U[[1],:] = [x0,y0]
    
    for i = 1:N-1
        t[i+1] = t[i] + h
        m1 = F(U[[i],:],t[i])
        m2 = F(U[[i],:] .+ transpose(h*m1), t[i])
        U[[i+1],:] = U[[i],:] + transpose(0.5*h*(m1 + m2))
    end
    return U,t
end
