from EvolAlgorithm import EvolAlgorithm
from Popul import PopulInt, PopulReal
from MotifFinding import MotifFinding
from MyMotifs import MyMotifs


def createMatZeros(nl, nc):
    res = []
    for _ in range(0, nl):
        res.append([0]*nc)
    return res


def printMat(mat):
    for i in range(0, len(mat)):
        for j in range(len(mat[i])):
            print(f"{mat[i][j]:.3f}", end=' ')
        print()


class EAMotifsInt (EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename, "dna")
        indsize = len(self.motifs)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulInt(self.popsize, indsize,
                              maxvalue, [])

    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            fit = self.motifs.score(sol)
            ind.setFitness(fit)


class EAMotifsReal(EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename, "dna")
        indsize = self.motifs.motifSize * len(self.motifs.alphabet)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulReal(self.popsize, indsize, maxvalue, [])

    def listToPwm(self, l):                #leva como parâmtero uma lista
        noSymb = len(self.motifs.alphabet) #nº de símbolos do alfabeto
        tamanhoMotif = self.motifs.motifSize    #tamanho do motif

        pwm = createMatZeros(noSymb, tamanhoMotif) #cria uma matriz de zeros com nº de linhas = nº de símbolos
                                              #e nº de colunas = tamanho do motif

        for i in range(0, len(l), noSymb):  #iterar sobre a lista l, com um intervalo = nº de símbolos
            index_coluna = i / noSymb
            col = l[i:i + noSymb]    #criar uma coluna com os elementos da lista com um intervalo = nº símbolos
            soma_coluna = sum(col)
            for j in range(noSymb): #por cada linha j
                pwm[j][int(index_coluna)] = col[j] / soma_coluna #estamos a preencher a matriz PWM por linhas

        # EXEMPLO:
        # [1,2,3,4,5,6,7,8]

        # A 1
        # C 2
        # T 3
        # G 4

        # A 1/10 2/10 3/10 4/10
        # C
        # T
        # G
    def evaluate(self, indivs):
        for i in range(len(indivs)):                    #vamos iterar sobre cada indivíduo
            ind = indivs[i]
            solucao = ind.getGenes()                    #obter a lista com os genes do indivíduo
            self.motifs.pwm = self.listToPwm(solucao)   #criar a matriz PWM com a lista de genes
            s = []                                      #lista com indices
            for seq in self.motifs.seqs:                #para cada sequência
                p = self.motifs.mostProbableSeq(seq)    #p vai ser o índice do motif mais provável na sequência
                s.append(p)                             #coloca-se o valor p na lista s
            fit = self.motifs.score(s)                 #calcula o score da lista s, o qual vai ser o valor fitness
            ind.setFitness(fit)                        #guarda o valor fitness para cada indivíduo

def test1():
    ea = EAMotifsInt(100, 1000, 50, "exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


def test2():
    ea = EAMotifsReal(100, 2000, 50, "exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


#test1()
test2()
