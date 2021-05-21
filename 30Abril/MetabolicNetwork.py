# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class MetabolicNetwork (MyGraph):#sub-classe da classe MyGraph (herança)-> deriva de uma classe
    
    def __init__(self, network_type = "metabolite-reaction", split_rev = False):
        #split_rev-> diz se representa reações reversiveis ou não, ou seja, se nos mostra as reações reversiveis atraves do grafo
        #indica se se consideram as reações reversíveis como duas reações distintas, uma para cada direção (True) ou não (False)
        #Tipo de rede-> é o que lhe passarmos
        MyGraph.__init__(self, {})#chamar o construtor da classe MyGraph (classe maior) -> estamos a chama-lo com o dicionario vazio para construir o grafo
        self.net_type = network_type #tipo de rede: “metabolite-reaction”, “reaction”, “metabolite”
        self.node_types = {}##guarda dicionário indicando listas de nós de cada tipo (no caso de ser uma rede “metabolite-reaction” – grafo bipartido)
        if network_type == "metabolite-reaction":#se for do tipo metabolito-reacao
            self.node_types["metabolite"] = []#estou a dizer que tenho nós deste tipo dentro da network
            self.node_types["reaction"] = []#estou a dizer que tenho nós deste tipo dentro da network
        self.split_rev =  split_rev
    
    def add_vertex_type(self, v, nodetype):#adicionar um nó/vertice
        self.add_vertex(v)#chama o add_vertex da classe mãe e cria um nó
        self.node_types[nodetype].append(v)#adiciona o no ao tipo de nó correspondente (metabolito ou reacao)
    
    def get_nodes_type(self, node_type):#obter o tipo de nós
        if node_type in self.node_types:#ver todos os nós no dicionario
            return self.node_types[node_type]#dar todos os values do nó, ou seja, todos os arcos/ligações entre aquele nó e outros nós
        else: return None
    
    def load_from_file(self, filename):#load do ficheiro para ter as reações e criar o grafo
        rf = open(filename)
        gmr = MetabolicNetwork("metabolite-reaction")#metabolito-reacao
        for line in rf:#ler cada linha do ficheiro
            if ":" in line:
                tokens = line.split(":")
                reac_id = tokens[0].strip()
                gmr.add_vertex_type(reac_id, "reaction")
                rline = tokens[1]
            else: raise Exception("Invalid line:")                
            if "<=>" in rline:
                left, right = rline.split("<=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_vertex_type(reac_id+"_b", "reaction")
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id+"_b", met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_edge(met_id, reac_id+"_b")
                        gmr.add_edge(reac_id, met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
            elif "=>" in line:
                left, right = rline.split("=>")#irreversivel
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(met_id, reac_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(reac_id, met_id)
            else: raise Exception("Invalid line:")    

        
        if self.net_type == "metabolite-reaction": #se for metabolito-reacao
            self.graph = gmr.graph#vai ser o grafo que está no objeto gmr em cima
            self.node_types = gmr.node_types
        elif self.net_type == "metabolite-metabolite":
            self.convert_metabolite_net(gmr)#teho que converter o grafo gmr em um grafo só de metabolitos
        elif self.net_type == "reaction-reaction": 
            self.convert_reaction_graph(gmr)#tenho de converter o grafo gmr em um grafo só de reações
        else: self.graph = {}
        
        
    def convert_metabolite_net(self, gmr):#criar a rede metabolito-metabolito
        for m in gmr.node_types["metabolite"]:#para cada metabolito em gmr
            self.add_vertex(m)#adicionar um novo metabolito ao novo grafo
            sucs = gmr.get_successors(m)#sucessores de m
            for s in sucs:#para cada sucessor de m
                sucs_r = gmr.get_successors(s)#vamos ver o sucessor dos sucessores de m:
                #por exemplo: M1 tem um sucessor R1, então vamos ver os sucessores de R1 que são M3 e M4
                for s2 in sucs_r:#para cada sucessor-sucessor
                    if m != s2:#se o m for diferente do sucessor-sucessor
                        self.add_edge(m, s2)#adicionar um arco de metabolitos

        
    def convert_reaction_graph(self, gmr):#criar a rede reação-reação
        for r in gmr.node_types["reaction"]:#para cada reação em gmr
            self.add_vertex(r)#adicionar uma nova reação ao novo grafo
            sucs = gmr.get_successors(r)#sucessores de r (reação)
            for s in sucs:#para cada sucessor de r
                sucs_r = gmr.get_successors(s)#vamos ver o sucessor dos sucessores de r:
                #por exemplo: R1 tem dois sucessores M3 e M4, então vamos ver os sucessores de M4 que são R2 e R3 e de M3 que não tem
                for s1 in sucs_r:#para cada sucessor-sucessor
                    if r != s1:#se o r não for igual ao sucessor-sucessor
                        self.add_edge(r, s1)#adicionar um arco de reações


def test1():
    m = MetabolicNetwork("metabolite-reaction")#cria o tipo de rede
    m.add_vertex_type("R1","reaction")#adicionar nó
    m.add_vertex_type("R2","reaction")
    m.add_vertex_type("R3","reaction")
    m.add_vertex_type("M1","metabolite")
    m.add_vertex_type("M2","metabolite")
    m.add_vertex_type("M3","metabolite")
    m.add_vertex_type("M4","metabolite")
    m.add_vertex_type("M5","metabolite")
    m.add_vertex_type("M6","metabolite")
    m.add_edge("M1","R1")#adicionar arco entre os nos (ao nó M1 vou adicionar uma ligação R1)
    m.add_edge("M2","R1")
    m.add_edge("R1","M3")
    m.add_edge("R1","M4")
    m.add_edge("M4","R2")
    m.add_edge("M6","R2")
    m.add_edge("R2","M3")
    m.add_edge("M4","R3")
    m.add_edge("M5","R3")
    m.add_edge("R3","M6")
    m.add_edge("R3","M4")
    m.add_edge("R3","M5")
    m.add_edge("M6","R3")
    m.print_graph()
    print("Reactions: ", m.get_nodes_type("reaction") )#dizer qual e a lista de reações
    print("Metabolites: ", m.get_nodes_type("metabolite") )#dizer qual é a lista de metebolitos

        
def test2():
    print("metabolite-reaction network:")
    mrn = MetabolicNetwork("metabolite-reaction")
    mrn.load_from_file("example-net.txt")
    mrn.print_graph()
    print("Reactions: ", mrn.get_nodes_type("reaction") )
    print("Metabolites: ", mrn.get_nodes_type("metabolite") )
    print()
    
    print("metabolite-metabolite network:")
    mmn = MetabolicNetwork("metabolite-metabolite")
    mmn.load_from_file("example-net.txt")
    mmn.print_graph()
    print()
    
    print("reaction-reaction network:")
    rrn = MetabolicNetwork("reaction-reaction")
    rrn.load_from_file("example-net.txt")
    rrn.print_graph()
    print()
    
    print("metabolite-reaction network (splitting reversible):")
    mrsn = MetabolicNetwork("metabolite-reaction", True)
    mrsn.load_from_file("example-net.txt")
    mrsn.print_graph()
    print()
    
    print("reaction-reaction network (splitting reversible):")
    rrsn = MetabolicNetwork("reaction-reaction", True)
    rrsn.load_from_file("example-net.txt")
    rrsn.print_graph()
    print()

def test3():
    mrn = MetabolicNetwork("metabolite-reaction")
    mrn.load_from_file("ecoli.txt")
    #mrn.print_graph()
    #print("Reactions: ", mrn.get_nodes_type("reaction") )
    #print("Metabolites: ", mrn.get_nodes_type("metabolite") )
    #print()
    print (mrn.mean_degree("out"))
    d = mrn.prob_degree("out")# dicionario
    for x in sorted(d.keys()):#para cada key no dicionario d por ordem crescente
        print (x, 't', d[x])#grau, tempo, value dessa key(grau)

#test1()
print()
#test2()
test3()
