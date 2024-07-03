include("matrices.jl")
using LinearAlgebra

#Quociente de rayleigh
function rayleigh(A,x)
    return (x'*A*x)/(x'*x);
end

#Função do método das potencias sem shifts para achar o maior autovalor em módulo
function potencias(A,x0=ones(size(A,1)),maxIt=10000,tol=1.0e-6)
    x=x0; x1=x; it=0
    while it<maxIt
        x=A*x; x=x/norm(x)
        #Critério de parada
        if norm(x1-x)<tol
            break
        end
        x1=x; it+=1
    end
    return x, rayleigh(A,x), it
end

#Função do método das potencias com shifts para achar o maior autovalor em módulo
function potenciasShift(A,x0=ones(size(A,1)), maxIt=10000,tol=1.0e-6)
    n=size(A,1); it=0
    #Um tanto de iterações do método das potencias usual para obter x0 para o método -> Com poucas iterações essa estimativa é uma merda
    x,ray,_=potencias(A,x0,30); x1=x
    while it<maxIt && det(A-ray*I)!=0
        x=(A-ray*I)\x; x=x/norm(x)
        #Critério de parada
        if norm(x1-x)<tol
            break
        end
        x1=x; ray=rayleigh(A,x); it+=1
    end
    return x, rayleigh(A,x), it
end

#Função do método de potenciais inverso para achar o menor autovalor em módulo
function potenciasInv(A,x0=ones(size(A,1)),maxIt=10000,tol=1.0e-6)
    x=x0; x1=x; it=0
    while it<maxIt
        x=A\x; x=x/norm(x)
        #Critério de parada
        if norm(x1-x)<tol
            break
        end
        x1=x; it+=1
    end
    return x, rayleigh(A,x), it
end

#Função do algoritmo QR
function autoQR(A,maxIt=10000,tol=1.0e-10)
    n=size(A,1); it=0
    v=Matrix{Float64}(I,n,n) #Matriz dos autovetores
    while it<maxIt
        Q, R=qr(A)
        A=R*Q; v*=Q
        #Critério de parada
        if maximum(abs.(tril(A,-1)))<tol
            break
        end
        it+=1
    end
    val=diag(A); #Autovalores
    return val, v
end