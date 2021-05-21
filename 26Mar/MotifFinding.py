# -*- coding: utf-8 -*-
"""

"""

from MySeq import MySeq
from MyMotifs import MyMotifs


class MotifFinding:

    def __init__(self, size=8, seqs=None):
        self.motifSize = size
        if (seqs != None):
            self.seqs = seqs
            self.alphabet = seqs[0].alfabeto()
        else:
            self.seqs = []

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, n):
        return self.seqs[n]

    def seqSize(self, i):
        return len(self.seqs[i])

    def readFile(self, fic, t):
        for s in open(fic, "r"):
            self.seqs.append(MySeq(s.strip().upper(), t))
        self.alphabet = self.seqs[0].alfabeto()

    def createMotifFromIndexes(self, indexes):
        pseqs = []
        for i, ind in enumerate(indexes):
            pseqs.append(MySeq(self.seqs[i][ind:(ind + self.motifSize)], self.seqs[i].tipo))
        return MyMotifs(pseqs)

    # SCORES

    def score(self, s):
        score = 0
        motif = self.createMotifFromIndexes(s)
        motif.doCounts()
        mat = motif.counts
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score += maxcol
        return score

    def score_pseudo(self, s):
        score = 0
        motif = self.createMotifFromIndexes(s)
        motif.doCounts_pseudo()    #só é preciso alterar aqui a função doCounts
        mat = motif.counts         #da classe MyMotifs para doCounts_pseudo
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score += maxcol
        return score

    def scoreMult(self, s):
        score = 1.0
        motif = self.createMotifFromIndexes(s)
        motif.createPWM()
        mat = motif.pwm
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score *= maxcol
        return score

        # EXHAUSTIVE SEARCH

    def nextSol(self, s):
        nextS = [0] * len(s)
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if (pos < 0):
            nextS = None
        else:
            for i in range(pos):
                nextS[i] = s[i]
            nextS[pos] = s[pos] + 1;
            for i in range(pos + 1, len(s)):
                nextS[i] = 0
        return nextS

    def exhaustiveSearch(self):
        melhorScore = -1
        res = []
        s = [0] * len(self.seqs)
        while (s != None):
            sc = self.score(s)
            if (sc > melhorScore):
                melhorScore = sc
                res = s
            s = self.nextSol(s)
        return res

    # BRANCH AND BOUND

    def nextVertex(self, s):
        res = []
        if len(s) < len(self.seqs):  # internal node -> down one level
            for i in range(len(s)):
                res.append(s[i])
            res.append(0)
        else:  # bypass
            pos = len(s) - 1
            while pos >= 0 and s[pos] == self.seqSize(pos) - self.motifSize:
                pos -= 1
            if pos < 0:
                res = None  # last solution
            else:
                for i in range(pos): res.append(s[i])
                res.append(s[pos] + 1)
        return res

    def bypass(self, s):
        res = []
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if pos < 0:
            res = None
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos] + 1)
        return res

    def branchAndBound(self):
        melhorScore = -1
        melhorMotif = None
        size = len(self.seqs)
        s = [0] * size
        while s != None:
            if len(s) < size:
                optimScore = self.score(s) + (size - len(s)) * self.motifSize
                if optimScore < melhorScore:
                    s = self.bypass(s)
                else:
                    s = self.nextVertex(s)
            else:
                sc = self.score(s)
                if sc > melhorScore:
                    melhorScore = sc
                    melhorMotif = s
                s = self.nextVertex(s)
        return melhorMotif

    # Consensus (heuristic)

    def heuristicConsensus(self):
        # procurar as posições para o motif nas 2 primeiras sequencias
        # procura exaustiva nas 2 primeiras sequencias
        mf = MotifFinding(self.motifSize, self.seqs[:2])
        s = mf.exhaustiveSearch()  # resultado, vamos supor, (1,3)
        # avaliar a melhor posição para cada uma das sequencias seguintes uma a uma,
        # guardando a melhor posição (maximiza o score)
        for i in range(2, len(self.seqs)):
            s.append(0)
            melhorScore = -1
            melhorPosicao = 0
            for j in range(self.seqSize(i) - self.motifSize + 1):
                s[i] = j
                score_atual = self.score(s)
                if score_atual > melhorScore:
                    melhorScore = score_atual
                    melhorPosicao = j
                s[i] = melhorPosicao
        return s

    # Consensus (heuristic)

    def heuristicStochastic(self):
        from random import randint
        s = [0] * len(self.seqs)
        # passo1: inicia todas as posiçoes com valores aleatorios
        for i in range(len(self.seqs)):
            s[i] = randint(0, self.seqSize(i) - self.motifSize)
        best_score = self.score(s)
        improve = True
        while improve:
            # passo2: criar o perfil
            # controi o perfil com base nas posiçoes iniciais s
            motif = self.createMotifFromIndexes(s)
            motif.createPWM()
            # passo 3
            # avalia a melhor posiçao inicial com base no perfil
            for i in range(len(self.seqs)):
                s[i] = motif.mostProbableSeq(self.seqs[i])
            # passo 4
            # verifica se houve melhoria ou não. se não houver pára o loop
            scr = self.score(s)
            if scr > best_score:
                best_score = scr
            else:
                improve = False
        return s
'''
heuristicStochastic com pseudo contagens
'''
    def heuristicStochastic_pseudo(self):
        from random import randint
        s = [0] * len(self.seqs)
        # passo1: inicia todas as posiçoes com valores aleatorios
        for i in range(len(self.seqs)):
            s[i] = randint(0, self.seqSize(i) - self.motifSize)

        best_score = self.score_pseudo(s)
        improve = True
        while improve:
            # passo2: criar o perfil
            # controi o perfil com base nas posiçoes iniciais s
            motif = self.createMotifFromIndexes(s)
            motif.createPWM_pseudo()
            # passo 3
            # avalia a melhor posiçao inicial com base no perfil
            for i in range(len(self.seqs)):
                s[i] = motif.mostProbableSeq(self.seqs[i])
            # passo 4
            # verifica se houve melhoria ou não. se não houver para o loop
            scr = self.score_pseudo(s)
            if scr > best_score:
                best_score = scr
            else:
                improve = False
        return s

    # Gibbs sampling

    def gibbs(self):
        from random import randint
        s = [0] * len(self.seqs)
        # passo1: inicia todas as posiçoes com valores aleatorios
        for i in range(len(self.seqs)):                         #s é uma lista com valores aleatórios para 
            s[i] = randint(0, self.seqSize(i) - self.motifSize) #as posições iniciais dos motifs nas sequências                                                         
        best_score = self.score(s) #soma dos maiores valores das células da matriz de contagens de cada coluna                            
        improve = True
        while improve:
            # passo 2
            # escolher uma das sequencias aleatoriamente
            seq_rmv = randint(0, len(self.seqs) - 1)
            # passo 3
            # removemos da lista de sequências a sequencia escolhida
            seq = self.seqs.pop(seq_rmv) #pop tira da lista a sequência através da posição desta na lista de seqs
                                         #e retorna a sequência removida na variável "seq"
            s_partial = s.copy()   #copiamos a lista s para a variável "s_partial"
            s_partial.pop(seq_rmv) #para depois removermos a posição inicial correspondente
                                   #à posição da sequência removida da lista "s_partial"
            motif = self.createMotifFromIndexes(s_partial) #motif vai ser um objecto da classe MyMotifs (MyMotifs(pseqs)) em que pseqs
                                                           #é uma lista com os motifs correspondentes às posições iniciais que estão em "s_partial"
            motif.createPWM() #criamos o perfil PWM (esta não vai ter a sequência que foi removida)
            # melhor posição inicial na seq considerando o perfil
            s[seq_rmv] = motif.mostProbableSeq(seq) #índice do segmento na sequência que deu uma probabilidade maior
            '''#p = motif.probAllPositions(seq)
            #position = self.roulette(p)
            #s.insert(seq_rmv, position)'''
            self.seqs.insert(seq_rmv, seq) #inserimos a sequência removida na lista de sequências na posição original desta
            scr = self.score(s) #calcula o novo score
            # verificar se houve melhorias
            if scr > best_score:
                best_score = scr #se o novo score for melhor que o inicial, então voltamos a iterar o bloco while
            else:
                improve = False #se não for melhor, significa que já temos o melhor score e não voltamos a correr o bloco while

        return s #retorna a lista de posições iniciais

    '''EA algorithms
    
    def gibbs(self, numits):
        pass
        
    '''
    def gibbs_pseudo(self):
        from random import randint
        s = [0] * len(self.seqs)
        # passo1: inicia todas as posiçoes com valores aleatorios
        for i in range(len(self.seqs)):                         #s é uma lista com valores aleatórios para 
            s[i] = randint(0, self.seqSize(i) - self.motifSize) #as posições iniciais dos motifs nas sequências                                                         
        best_score = self.score_pseudo(s) #soma dos maiores valores das células da matriz de contagens de cada coluna                            
        improve = True
        while improve:
            # passo 2
            # escolher uma das sequencias aleatoriamente
            seq_rmv = randint(0, len(self.seqs) - 1)
            # passo 3
            # removemos da lista de sequências a sequencia escolhida
            seq = self.seqs.pop(seq_rmv) #pop tira da lista a sequência através da posição desta na lista de seqs
                                         #e retorna a sequência removida na variável "seq"
            s_partial = s.copy()   #copiamos a lista s para a variável "s_partial"
            s_partial.pop(seq_rmv) #para depois removermos a posição inicial correspondente
                                   #à posição da sequência removida da lista "s_partial"
            motif = self.createMotifFromIndexes(s_partial) #motif vai ser um objecto da classe MyMotifs (MyMotifs(pseqs)) em que pseqs
                                                           #é uma lista com os motifs correspondentes às posições iniciais que estão em "s_partial"
            motif.createPWM_pseudo() #criamos o perfil PWM (esta não vai ter a sequência que foi removida)
            # melhor posição inicial na seq considerando o perfil
            s[seq_rmv] = motif.mostProbableSeq(seq) #índice do segmento na sequência que deu uma probabilidade maior
            '''#p = motif.probAllPositions(seq)
            #position = self.roulette(p)
            #s.insert(seq_rmv, position)'''
            self.seqs.insert(seq_rmv, seq) #inserimos a sequência removida na lista de sequências na posição original desta
            scr = self.score_pseudo(s) #calcula o novo score
            # verificar se houve melhorias
            if scr > best_score:
                best_score = scr #se o novo score for melhor que o inicial, então voltamos a iterar o bloco while
            else:
                improve = False #se não for melhor, significa que já temos o melhor score e não voltamos a correr o bloco while

        return s #retorna a lista de posições iniciais 

    def roulette(self, f):
        from random import random
        tot = 0.0
        for x in f: tot += (0.01 + x)
        val = random() * tot
        acum = 0.0
        ind = 0
        while acum < val:
            acum += (f[ind] + 0.01)
            ind += 1
        return ind - 1


# tests

def test1():
    sm = MotifFinding()
    sm.readFile("exemploMotifs.txt", "dna")
    sol = [25, 20, 2, 55, 59]
    sa = sm.score(sol)
    print(sa)
    scm = sm.scoreMult(sol)
    print(scm)


def test2():
    print("Test exhaustive:")
    seq1 = MySeq("ATAGAGCTGA", "dna")
    seq2 = MySeq("ACGTAGATGA", "dna")
    seq3 = MySeq("AAGATAGGGG", "dna")
    mf = MotifFinding(3, [seq1, seq2, seq3])
    sol = mf.exhaustiveSearch()
    print("Solution", sol)
    print("Score: ", mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

    print("Branch and Bound:")
    sol2 = mf.branchAndBound()
    print("Solution: ", sol2)
    print("Score:", mf.score(sol2))
    print("Consensus:", mf.createMotifFromIndexes(sol2).consensus())

    print("Heuristic consensus: ")
    sol1 = mf.heuristicConsensus()
    print("Solution: ", sol1)
    print("Score:", mf.score(sol1))


def test3():
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt", "dna")
    print("Branch and Bound:")
    sol = mf.branchAndBound()
    print("Solution: ", sol)
    print("Score:", mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())


def test4():
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt", "dna")
    print("Heuristic stochastic")
    sol = mf.heuristicStochastic()
    print("Solution: ", sol)
    print("Score:", mf.score(sol))
    print("Score mult:", mf.scoreMult(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

    print("Gibbs")
    sol2 = mf.gibbs()
    print("Score:", mf.score(sol2))
    print("Score mult:", mf.scoreMult(sol2))
    print("Score mult:", mf.scoreMult(sol2))
    print("Consensus:", mf.createMotifFromIndexes(sol2).consensus())

def test5():
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt", "dna")
    print("PSEUDO Heuristic stochastic")
    sol = mf.heuristicStochastic_pseudo()
    print("Solution: ", sol)
    print("Score:", mf.score(sol))
    print("Score mult:", mf.scoreMult(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

    print("PSEUDO Gibbs")
    sol2 = mf.gibbs_pseudo()
    print("Score:", mf.score(sol2))
    print("Score mult:", mf.scoreMult(sol2))
    print("Score mult:", mf.scoreMult(sol2))
    print("Consensus:", mf.createMotifFromIndexes(sol2).consensus())


if __name__ == '__main__':
    test1()
    test2()
    test3()
    test4()
    test5()
