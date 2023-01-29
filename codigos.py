import numpy as np

def trocaLinhaB(b,inicio,final):
    listaAux = b[inicio].copy()
    b[inicio] = b[final]
    b[final] = listaAux
    
def trocaColunaB(b,inicio,final):
    listaAux = b[:,inicio].copy()
    b[:,inicio] = b[:,final]
    b[:,final] = listaAux
    
def fazPivoteamento(M, m, n): # matriz M mxn
    primeira = 1;
    linha = 0;
    for j in range(n): # caminha nas colunas esquerda para a direita
        for i in range(m): # caminha nas linhas cima para baixo
            if(M[i,j] == 1): # procura valor 1
                if(primeira == 1 and i >= linha): # se for o primeiro valor 1 da coluna entra
                    trocaLinhaB(M,i,linha)   # coloca a linha i para a posição N da matriz identidade
                    primeira = 0;
                    aux = (linha-1)

                    while aux>=0:   # irá verificar as linhas acima da linha N
                        if(M[aux,j] == 1):
                            M[aux,:] = (M[aux,:] + M[linha,:])%2
                        aux = aux-1
                    linha=linha+1;

                elif(primeira==0): # se NÃO for o primeiro valor 1 da coluna entra
                    M[i,:] = (M[i,:] + M[linha-1,:])%2
        primeira = 1;
  return M;
  
def ordenaMatrizB(B):
    M = B.copy()

    nlinhas = len(M)
    ncolunas = len(M[0])
    columns = np.arange(1,ncolunas+1)

    posto = nlinhas;
    diag = 0; 

    b=np.cumsum(M,axis=1)

    for i in range(nlinhas):
        if (M[i,:].sum()==0):
            posto = posto-1
    while diag < posto :
        if (M[:,diag].sum() == 1 ):
            diag = diag+1
        else:
            aux = diag
            while (aux < posto):
                aux = aux+1
                if (M[:,aux].sum() == 1):
                    trocaColunaB(M,diag,aux)
                    trocaLinhaB(columns,diag,aux)
                    diag = diag+1
  return M, columns, posto;
  
  def Rotina(E):
    nlinhas = len(E)
    ncolunas = len(E[0])

    # Passo 1, recebendo recebendo a matriz de ciclo e fazendo o pivoteamento

    MatrizPivoteada = fazPivoteamento(E, nlinhas, ncolunas)
    C, variaveis, posto = ordenaMatrizB(MatrizPivoteada)

    #Este bloco controla os indices da matriz de ciclo fundamental
    IF = variaveis[:posto]
    CF = variaveis[posto:len(C[0])]
    variaveis = np.concatenate((CF,IF), axis=0)

    # Este bloco irá remover a matriz identidade e fazer a transposta de Cf e juntar com a identidade e encontrar B
    i = 0
    while i<posto:
      C = np.delete(C, 0, axis = 1)
      i+=1
    while (len(C) > posto):
      C = np.delete(C, len(C)-1, axis = 0)
    B = C.transpose()

    IB=np.identity(len(B), dtype=int)
    C = np.concatenate((IB,B), axis=1)

    # Passo 3, a partir da matriz de corte fundamental, encontrar a matriz de incidência
    B = C.copy()
    i=0

    # Este bloco fará operações para que não tenham colunas com mais de 2 valores não nulos
    while i<len(B[0]):
        if (B[:,i].sum()>2):
          vetor = B[:,i]
          primeiro = 1
          acabou = 0
          a=0
          while a < len(vetor) and acabou == 0:
            fixo = 0
            if vetor[a] == 1 and primeiro == 1:
              fixo = a
              primeiro = 0
            elif vetor[a] == 1:
              B[a,:] = (B[fixo,:] + B[a,:])%2
              vetor = B[:,i]
              a = a + 1 
            elif vetor.sum() == 2:
              acabou = 1
            a = a + 1 
          i=i+1
        else:
          i=i+1

    # Este bloco irá adicionar a ultima linha de Q
    novoVetor = np.zeros(len(B[0]), dtype=int)
    for i in range(len(B[0])):
      if(B[:,i].sum()<2):
        novoVetor[i] = 1
    aux = np.array([novoVetor], dtype=int)
    Q = np.concatenate((B,aux))
    return Q, variaveis
