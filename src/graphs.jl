include("metodos.jl")
include("matrices.jl")
using LinearAlgebra, DataFrames, PrettyTables, Suppressor, Printf, Statistics

#Função para rodar os testes, devolve os erros relativos e tempo de cada metodo
function run(dims)
    n=length(dims)
    #Alocando variaveis para armazenar erro e tempo para os métodos das potencias usual
    errosMaior=zeros(n,3,3); temposMaior=zeros(n,3,3); valMax=zeros(3)
    #Alocando para o método das potencias inverso, para achar o menor
    errosMenor=zeros(n,2,3); temposMenor=zeros(n,2,3); valMin=zeros(2)

    #Rodando para cada caso
    for opt=1:3
        for i=1:n
            #Autovalores exatos
            autoval=reverse(diag(Mdiag(opt,dims[i])))

            #Alocando a matriz
            A=testMatrix(opt,dims[i]) 

            #Soluçao usando cada método das potencias direto
            temposMaior[i,1,opt]=@elapsed _, valMax[1]=potencias(A)
            temposMaior[i,2,opt]=@elapsed _, valMax[2]=potenciasShift(A)
            #Pelo método QR
            temposMaior[i,3,opt]=@elapsed _, qr_autoval=autoQR(A);
            temposMenor[i,2,opt]=temposMaior[i,3,opt]
            valMin[2]=qr_autoval[dims[i]]; valMax[3]=qr_autoval[1];
            #E pelo método das potencias inverso
            temposMenor[i,1,opt]=@elapsed _, valMin[1]=potenciasInv(A)

            #preenchendo os erros para o maior autovalor
            for j=1:3
                errosMaior[i,j,opt]=abs(valMax[j]-autoval[1])/abs(autoval[1])
            end

            #Preenchendo os erros para o menor autovalor
            for j=1:2
                errosMenor[i,j,opt]=abs(valMin[j]-autoval[dims[i]])/abs(autoval[dims[i]])
            end
        end
    end

    return errosMaior, temposMaior, errosMenor, temposMenor
end

#Função para montar as tabelas
function tables_autoval(dim)
    eM, tM, em, tm=run(dim) #Executando o teste 

    #Formatando as tabelas para dataframes
    dfeM=[hcat(dim,eM[:,:,i]) for i=1:3] ; dftM=[hcat(dim,tM[:,:,i]) for i=1:3] 
    dfem=[hcat(dim,em[:,:,i]) for i=1:3]; dftm=[hcat(dim,tm[:,:,i]) for i=1:3]

    #Passando cada matriz para dataframe
    dt_dfeM=[]; dt_dftM=[]; dt_dfem=[]; dt_dftm=[]
    header=["n", "Pot.", "Pot. Shifts", "QR"]
    for i=1:3
        push!(dt_dfeM, DataFrame(dfeM[i],Symbol.(header)))
        push!(dt_dftM, DataFrame(dftM[i],Symbol.(header)))
        #Passando a coluna do N para inteiro 
        dt_dfeM[i].n=Int.(dt_dfeM[i].n); dt_dftM[i].n=Int.(dt_dftM[i].n)
        #Passando cada entrada nas células para notação cientifica
        for col in names(dt_dfeM[i])[2:end]
            dt_dfeM[i][!, col]=string.(@sprintf("%.3e", x) for x in dt_dfeM[i][!,col])
            dt_dftM[i][!, col]=string.(@sprintf("%.3e", x) for x in dt_dftM[i][!,col])
        end
    end
    header=["n", "Pot. Inv", "QR"]
    for i=1:3
        push!(dt_dfem, DataFrame(dfem[i],Symbol.(header)))
        push!(dt_dftm, DataFrame(dftm[i],Symbol.(header)))
        #Passando a coluna do N para inteiro
        dt_dfem[i].n=Int.(dt_dfem[i].n); dt_dftm[i].n=Int.(dt_dftm[i].n)
        #Passando cada entrada nas células para notação cientifica
        for col in names(dt_dfem[i])[2:end]
            dt_dfem[i][!, col]=string.(@sprintf("%.3e", x) for x in dt_dfem[i][!,col])
            dt_dftm[i][!, col]=string.(@sprintf("%.3e", x) for x in dt_dftm[i][!,col])
        end
    end

    #Printando e colocando em um arquivo de texto.
    open("output.txt","w") do io
        for i=1:3
            #Escrevendo os resultados para encontrar o maior autovalor 
            write(io,"Problema ("*Char(64+i)*") - Erros relativos - Maior Autovalor\n")
            latex_tab=@capture_out pretty_table(dt_dfeM[i],backend=Val(:latex),header=names(dt_dfeM[i]))
            write(io,latex_tab); write(io,'\n')
            write(io,"Problema ("*Char(64+i)*") - Tempos - Maior Autovalor\n");
            latex_tab=@capture_out pretty_table(dt_dftM[i],backend=Val(:latex),header=names(dt_dftM[i]))
            write(io,latex_tab); write(io,'\n')

            #E para o menor autovalor
            write(io,"Problema ("*Char(64+i)*") - Erros relativos - Menor Autovalor\n")
            latex_tab=@capture_out pretty_table(dt_dfem[i],backend=Val(:latex),header=names(dt_dfem[i]))
            write(io,latex_tab); write(io,'\n')
            write(io,"Problema ("*Char(64+i)*") - Tempos - Menor Autovalor\n")
            latex_tab=@capture_out pretty_table(dt_dftm[i],backend=Val(:latex),header=names(dt_dftm[i]))
            write(io,latex_tab); write(io,'\n')
        end
    end
end

#Função que devolve um dataframe com os erros obtidos na estimativa dos autovalores pelo algoritmo QR
function erroQR(opt,dim,nt)
    resultGeral=[]; result=[];
    for i=1:nt
        for n in dim
            A=testMatrix(opt,n) #Matriz de teste
            _, val=autoQR(A) #Autovalores obtidos pelo método
            autoval=reverse(diag(Mdiag(opt,n))) #Autovalores exatos
            push!(result,norm((autoval-val)/norm(autoval)))
        end
        push!(resultGeral,result)
    end
    resultGeral=[mean(v) for v in resultGeral]
    return mean(resultGeral)
end

dim=[2,5,10,100,1000]; tables_autoval(dim)