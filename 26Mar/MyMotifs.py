# -*- coding: utf-8 -*-
"""
@author: miguelrocha
"""

def createMatZeros (nl, nc):
    res = [ ] 
    for i in range(0, nl):
        res.append([0]*nc)
    return res

def printMat(mat):
    for i in range(0, len(mat)): print(mat[i])

class MyMotifs:

    def __init__(self, seqs):
        self.size = len(seqs[0])
        self.seqs = seqs # objetos classe MySeq
        self.alphabet = seqs[0].alfabeto()
        self.doCounts()
        self.createPWM()
        
    def __len__ (self):
        return self.size        
        
    def doCounts(self):
        self.counts = createMatZeros(len(self.alphabet), self.size)
        for s in self.seqs:
            for i in range(self.size):
                lin = self.alphabet.index(s[i])
                self.counts[lin][i] += 1
                
####doCounts, ou seja, criar matriz de contagens, mas com pseudocontagens
    def doCounts_pseudo(self):
        self.counts = createMatZeros(len(self.alphabet), self.size)
        for s in self.seqs:
            for i in range(self.size):
                lin = self.alphabet.index(s[i])
                self.counts[lin][i] += 1
        for i in range(len(self.counts)):           #len(self.counts)->no. de linhas
            for j in range(len(self.counts[0])):     #len(self.counts[0])->no. de colunas
                self.counts[i][j]+=1     #vamos percorrer as células por linha
                                                    #e adicionamos uma unidade ao valor de cada célula        
                
    def createPWM(self):
        if self.counts == None: self.doCounts()
        self.pwm = createMatZeros(len(self.alphabet), self.size)
        for i in range(len(self.alphabet)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j]) / len(self.seqs)

    def createPWM_pseudo(self):                                    #como alterámos os valores da matriz
        if self.counts == None: self.doCounts_pseudo()             #tivemos que criar um perfil probabilístico (position weighted matrix – PWM),    
        columnSum = 0                                              #em que, em vez de os valores das células serem iguais ao resultado da divisão
        for i in range(len(self.counts)):                          #dos valores das células da matriz de contagens (doCounts) pelo número total de sequências,
            columnSum += self.counts[i][0]                         #os valores das células da matriz de contagens vai ser dividido pela soma dos valores
        self.pwm = createMatZeros(len(self.alphabet), self.size)   #das colunas (columnSum) (coluna 0 -> soma-se os valores das células de todas as linhas, etc.) 
        for i in range(len(self.alphabet)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j]) / columnSum 
                
    def consensus(self):
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet) ):
                if self.counts[i][j] > maxcol: 
                    maxcol = self.counts[i][j]
                    maxcoli = i
            res += self.alphabet[maxcoli]        
        return res

    def maskedConsensus(self):
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet) ):
                if self.counts[i][j] > maxcol: 
                    maxcol = self.counts[i][j]
                    maxcoli = i
            if maxcol > len(self.seqs) / 2:
                res += self.alphabet[maxcoli]        
            else:
                res += "-"
        return res

    def probabSeq (self, seq):
        res = 1.0
        for i in range(self.size):
            lin = self.alphabet.index(seq[i])
            res *= self.pwm[lin][i]
        return res
    
    def probAllPositions(self, seq):
        res = []
        for k in range(len(seq)-self.size+1):
            res.append(self.probabSeq(seq))
        return res

    def mostProbableSeq(self, seq):
        maximo = -1.0
        maxind = -1
        for k in range(len(seq)-self.size):
            p = self.probabSeq(seq[k:k+ self.size])
            if(p > maximo):
                maximo = p
                maxind = k
        return maxind

def test():
    # test
    from MySeq import MySeq
    seq1 = MySeq("AAAGTT")
    seq2 = MySeq("CACGTG")
    seq3 = MySeq("TTGGGT")
    seq4 = MySeq("GACCGT")
    seq5 = MySeq("AACCAT")
    seq6 = MySeq("AACCCT")
    seq7 = MySeq("AAACCT")
    seq8 = MySeq("GAACCT")
    lseqs = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8]
    motifs = MyMotifs(lseqs)
    #printMat (motifs.counts)
    #printMat (motifs.pwm)
    #print(motifs.alphabet)
    
    print(motifs.probabSeq("AAACCT"))
    #print(motifs.probabSeq("ATACAG"))
    #print(motifs.mostProbableSeq("CTATAAACCTTACATC"))
    
    #print(motifs.consensus())
    #print(motifs.maskedConsensus())

if __name__ == '__main__':
    test()
