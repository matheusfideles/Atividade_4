include("metodos.jl")
include("matrices.jl")
using LinearAlgebra, Plots, PlotlyJS, DataFrames

#Função para rodar os testes, devolve os erros relativos e tempo de cada metodo
function run(dims)
    n=length(dims)
    #Alocando variaveis para armazenar erro e tempo para os métodos das potencias usual
    errosMaior=zeros(n,3,3); temposMaior=zeros(n,3,3); valMax=zeros(3)
    #Alocando para o método das potencias inverso, para achar o menor
    errosMenor=zeros(n,3,2); temposMenor=zeros(n,3,2); valMin=zeros(2)

    #Rodando para cada caso
    for opt=1:3
        for i=1:n
            #Autovalores exatos
            autoval=reverse(diag(Mdiag(opt,dims[i])))

            #Alocando a matriz
            A=testMatrix(opt,dims[i]) 

            #Soluçao usando cada método das potencias direto
            temposMaior[i,opt,1]=@elapsed _, valMax[1]=potencias(A)
            temposMaior[i,opt,2]=@elapsed _, valMax[2]=potenciasShift(A)
            #Pelo método QR
            temposMaior[i,opt,3]=@elapsed _, qr_autoval=autoQR(A); temposMenor[i,opt,2]=temposMaior[i,opt,3]
            valMin[2]=qr_autoval[dims[i]]; valMax[3]=qr_autoval[1];
            #E pelo método das potencias inverso
            temposMenor[i,opt,1]=@elapsed _, valMin[1]=potenciasInv(A)

            #preenchendo os erros para o maior autovalor
            for j=1:3
                errosMaior[i,opt,j]=abs(valMax[j]-autoval[1])/abs(autoval[1])
            end

            #Preenchendo os erros para o menor autovalor
            for j=1:2
                errosMenor[i,opt,j]=abs(valMin[j]-autoval[dims[i]])/abs(autoval[dims[i]])
            end
        end
    end

    return errosMaior, temposMaior, errosMenor, temposMenor
end

#Função para montar as tabelas
function tables_autoval(dim)
    eM, tM, em, tm=run(dim) #Executando o teste 

    #Formatando as tabelas para dataframes
    dfeM=[hcat(dim,eM[:,:,i]) for i=1:3] 
    dftM=[hcat(dim,tM[:,:,i]) for i=1:3] 
    dfem=[hcat(dim,em[:,:,i]) for i=1:2]
    dftm=[hcat(dim,tm[:,:,i]) for i=1:2]

    #Passando para latex e salvando
    

    #Para o met. das potencias inverso
    return dfeM,dftM,dfem,dftm
end

dim=[2,5,10,100,1000]