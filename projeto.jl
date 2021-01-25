using LinearAlgebra
function HouseholderQR(A)
    A = Float64.(A)
    k = 1
    b = 0
    u = [] ##u é o vetor que a gente vai construir para determinar Q
    auxiliar = 0
    m,n = size(A) ##m são as colunas e n as linhas de A
    if m != n 
        error("a matriz A deve ser nxn") 
    end
    while k < n
        u = [] ##vai resetar u quando terminar o loop
    ##cria o vetor u igual ao vetor v, onde v é a coluna de A
        for i in 1:m
            for j in k:k
                u = [u; A[i,j]]
            end
        end
        for i in 1:m
            if i < k
                u[i] = 0   ##todos os elementos de u acima do elemento k serão iguais a zero             
            end
        end
        for l in (k):m
            auxiliar += A[l,k]^2 ##auxiliar que calcula os elementos dentro da raíz da norma
        end
        s = sqrt(auxiliar) ##norma
        if s*A[k,k] <= 0 ##caso o produto de s e o elemento k sejam maiores que zero, ele muda o sinal de s
            s = -s
        end
        u[k] = A[k,k] - s    ## o elemento k do vetor u é igual o elemento k,k da matriz A
        b = 2/ (u'*u)  ##vai simplificar pra fazer a matriz Q
        k += 1
        Q = I - b*u*u' ##Define a matriz Q
        R = Q'*A*Q' ##Define a matriz R 
        println("\nA matriz Q$(k-1) = $Q'") ##Vai printar a matriz Q achada
        if Q'*R*Q == A  ##Caso QR = A, quebra o loop e vai para os prints.
            break
        end
    end
    Q = I - b*u*u'
    R = (Q'*A)
    #println("\nA matriz A = $A")
    println("\nMatriz R = $R")
    println("\nMatriz QR = $(Q'*R)")
    println("\nMatriz Q = $Q")
end

function GivensQR(a::Matrix)
    (m,n) = size(a)
    if m != n
        error("A matriz deve ser quadrada")
    end
    g = Matrix{Float64}(I,m,m)
    o = 1
    b = a
    u = 0
    t = 0
    p = 0
    if m == 2 
        g[i,j] = -b[i,j]/sqrt((b[i,j]^2)+(b[i-1,j]^2))
        g[i-1,j] = b[i-1,j]/sqrt((b[i-1,j]^2)+(b[i,j]^2))
        g[i-1,j+1] = b[i,j]/sqrt((b[i,j]^2)+(b[i-1,j]^2))
        g[i,j+1] = b[i-1,j]/sqrt((b[i-1,j]^2)+(b[i,j]^2))
        b = g*b
    else
        while o < ((m*m)-m) 
            for j=1:m
                for i=2:m #####se fosse 1:m ele nao conseguiria acessar os elementos
                    if i <= j #só vai acessar abaixo da diagonal principal 
                        u = u+1
                    elseif j==1 && i==m ############ elemento canto inferior esquerdo am1
                        g[i,j] = -(b[i,j]/sqrt((b[i,j]^2)+(b[1,1])^2))
                        g[i-2,j] = b[1,1]/sqrt((b[i,j]^2)+(b[1,1]^2))
                        g[i-1,j+1] = 1
                        g[i,j+2] = b[1,1]/sqrt((b[i,j]^2)+(b[1,1]^2))
                        g[i-2,j+2] = (b[i,j]/sqrt((b[i,j]^2)+(b[1,1])^2))
                        b = g*b
                        g = Matrix{Float64}(I,m,m)
                        t = t+1
                    else 
                        g[i,j] = -b[i,j]/sqrt((b[i,j]^2)+(b[i-1,j]^2))
                        g[i-1,j] = b[i-1,j]/sqrt((b[i-1,j]^2)+(b[i,j]^2))
                        g[i-1,j+1] = b[i,j]/sqrt((b[i,j]^2)+(b[i-1,j]^2))
                        g[i,j+1] = b[i-1,j]/sqrt((b[i-1,j]^2)+(b[i,j]^2))
                        b = g*b
                        g = Matrix{Float64}(I,m,m)
                        p = p+1
                    end
                    o = o+1
                end
            end
        end
    end
    l = b*inv(a)
    k = transpose(l)
    println("\nMatriz R = $b")
    println("\nMatriz Q = $k")
end