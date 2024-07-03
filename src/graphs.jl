include("metodos.jl")
include("matrices.jl")
using LinearAlgebra

#Função para plotar o gráfico de dispersao
function scatter_autoval()

end

#Função para rodar os testes, devolve os erros relativos e tempo de cada metodo
function run(dims)
    #Alocando variaveis para armazenar erro e tempo para os métodos das potencias usual
    errosMaior=zeros(dims,3,3); temposMaior=zeros(dims,3,3); tMax=zeros(3); valMax=zeros(3)
    #Alocando para o método das potencias inverso, para achar o menor
    errosMenor=zeros(dims,3,2); temposMenor=zeros(dims,3,2); tMin=zeros(2); valMin=zeros(2)

    #Rodando para cada caso
    for opt=1:3
        for i=1:length(dims)
            #Autovalores exatos
            autoval=Diag(Mdiag(opt,dims[i]))
            #Alocando a matriz
            A=testMatrix(opt,dims[i])
            #Soluçao usando cada método das potencias
            t[1]=@elapsed _, valMax[1]=potencias(A)
            t[2]=@elapsed _, valMax[2]=potenciasShift(A)
            t[3]=@elapsed _, valMax[3]=potenciasInv(A)
            #E pelo método QR
            t[4]=@elapsed _, qr_autoval=autoQR(A)
            ref_autoval=zeros(2);
            ref_autoval[1]=qr_autoval[1]; ref_autoval[2]=qr_autoval[dims[i]]
            #preenchendo os erros entre o maior autovalor pelo met. das pot e o exato
            for j=1:3
                errosMaior[i,opt,j]=[abs(val[j]-autoval[1])/abs(autoval[1])]
            end

            #Preenchendo os erros entre o menor autovalor pelo met. das pot e o exato
            erros[i,opt,4]=[]

        end
    end
end