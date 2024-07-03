using LinearAlgebra

#Retorna matriz diagonal da questão 1
function Mdiag(opt,n)
    D=zeros(n,n)
    if opt==1
        for i=1:n
            D[i,i]=(1.001)^i
        end
    elseif opt==2
        for i=1:n
            D[i,i]=(-1.001)^i
        end
    else
        for i=1:n
            D[i,i]=1.01^i
        end
    end
    return D
end

#Retorna matrix aleatoria com entradas entre -M e M
function randMatrix(n)
    M=999; R=M*(2*rand(n,n).-1)
    return R
end

#Constroi a matriz de teste
function testMatrix(opt,n)
    M=randMatrix(n); Q,_=qr(M)
    A=Q*Mdiag(opt,n)*Q'; A=Symmetric(A)
    return A
end 