# u_xx + u_yy - 10u = 0
# u = 0 at y = 1 and -1 <= x <= 1
# u = 1 at y = -1 -1 <= x <= 1
# u_x = -0.5u at x = 1 and -1 <= y <= 1
# u_x = 0.5u at -1 and -1 <= y <= 1

function helmholtz(m)
    h = 2/m
    d_up = ones(m-1)
    d_up[m-1] = 2
    
    d = (-4-10*h*h)*ones(m)
    d[1] = (-4-h-10*h*h)# du/dx = 0.5u
    d[m] = (-4+h-10*h*h)# du/dx = -0.5u
    
    d_down = ones(m-1)
    d_down[1] = 2
    
    A = Tridiagonal(d_down, d, d_up)
    A = kron(sparse(I,m,m),A)
    
    j=1
    for i=m+1:m*m
        A[i,j] = 1
        A[j,i] = 1
        j+=1
    end
    
    A = sparse(A)
    F = zeros(m*m)   
    F[1:m] .= -1
    
    return A\F,m,h,A
end
