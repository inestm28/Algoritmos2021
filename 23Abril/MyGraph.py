class MyGraph:

    def __init__(self, g={}):
        ''' Constructor - takes dictionary to fill the graph as input; default is empty dictionary '''
        self.graph = g

    def print_graph(self):
        ''' Prints the content of the graph as adjacency list '''
        for v in self.graph.keys():
            print(v, " -> ", self.graph[v])

    ## get basic info

    def get_nodes(self):
        ''' Returns list of nodes in the graph '''
        return list(self.graph.keys())

    def get_edges(self):
        ''' Returns edges in the graph as a list of tuples (origin, destination, peso) '''
        edges = []
        for v in self.graph.keys():
            for des in self.graph[v]:
                d, wg = des
                edges.append((v, d, wg))
        return edges

    def size(self):
        ''' Returns size of the graph : number of nodes, number of edges '''
        return len(self.get_nodes()), len(self.get_edges())

    ## add nodes and edges

    def add_vertex(self, v):
        ''' Add a vertex to the graph; tests if vertex exists not adding if it does '''
        if v not in self.graph.keys():
            self.graph[v] = []

    def add_edge(self, o, d, wg):
        ''' Add edge to the graph; if vertices do not exist, they are added to the graph '''
        if o not in self.graph.keys():  # verifica se ha o vertice o senao adiciona ao dicionario
            self.add_vertex(o)
        if d not in self.graph.keys():  # verifica se ha d vertice o senao adiciona ao dicionario
            self.add_vertex(d)
        des = []
        for j in self.graph[o]:
            d, iwg = j
            des.append(d)
        if d not in des:
            # verifica se ha ligação entre os dois vertices, caso contrario adiciona o vertice d à lista do vertice o
            self.graph[o].append((d, wg))

    ## successors, predecessors, adjacent nodes

    def get_successors(self, v):
        res = []
        for j in self.graph[v]:
            d, wg = j
            res.append(d)
        return res  # needed to avoid list being overwritten of result of the function is used

    def get_predecessors(self, v):
        lista = []
        for i in self.graph.keys():
            for tupls in self.graph[i]:
                d, wg = tupls
                if d == v:
                    lista.append(i)
        return lista

    def get_adjacents(self, v):
        suc = self.get_successors(v)  # buscar os sucessores
        pred = self.get_predecessors(v)  # buscar os predecessores
        res = pred
        for p in suc:  # adcionar os sucessores não presentes na lista
            if p not in res:
                res.append(p)
        return res

    ## degrees

    def out_degree(self, v):
        li = self.get_successors(v)
        return len(li)

    def in_degree(self, v):
        li = self.get_predecessors(v)
        return len(li)

    def degree(self, v):
        li = self.get_adjacents(v)
        return len(li)

    ## BFS and DFS searches

    def reachable_bfs(self, v):
        l = [v]  # sitio onde tou a começar, lista de coisas a começar, ou seja começa pelo no de origem
        res = []  # a lista do resultado de nos atingiveis
        while len(l) > 0:  # enquanto ha elementos na lista l (lista de queue)
            node = l.pop(0)  # isolar o 1º no na queue
            if node != v:  # se o v não entra na lista de nos acessiveis e for diferente é adiciona a lista
                res.append(node)  # controi o resultado sempre a colocar no fim, ou seja os primeiros ficam para ultimo
            for elem in self.graph[node]:  # correr os nos no node a pesquisar
                nwnode, wg = elem
                if nwnode not in res and nwnode not in l and nwnode != node:  # adicionar á queue
                    l.append(nwnode)
        return res

    def reachable_dfs(self, v):
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            s = 0
            if node != v:  # se o v não entra na lista de nos acessiveis e for diferente é adiciona a lista
                res.append(node)  # controi o resultado sempre a colocar no fim, ou seja os primeiros ficam para ultimo
            for elem in self.graph[node]:
                nwnode, wg = elem
                if nwnode not in res and nwnode not in l:
                    l.insert(s, nwnode)
                    s += 1
        return res

    def distance(self, s, d):
        if s == d:
            return 0
        else:
            l = [(s, 0)]  # sitio onde tou a começar, lista de coisas a começar, ou seja começa pelo no de origem
            vis = [s]  # a lista do resultado de nos atingiveis
            while len(l) > 0:  # enquanto ha elementos na lista l (lista de queue)
                node, swg = l.pop(0)  # isolar o 1º no na queue smw(sumaded weigh)
                for elem in self.graph[node]:  # correr os nos no node a pesquisar
                    nwnode, wg = elem
                    if nwnode == d: return swg + wg
                    if nwnode not in vis and nwnode not in l and nwnode != node:  # adicionar á queue
                        l.append((nwnode, swg + wg))
                        vis.append(nwnode)
            return None

    def shortest_path(self, s, d):  #algoritmo de djiskra
        if s == d:
            return [s, d]
        else:
            l = [(s, [], 0)]  # sitio onde tou a começar, lista de coisas a começar, ou seja começa pelo no de origem
            vis = [s]  # a lista do resultado de nos atingiveis
            while len(l) > 0:  # enquanto ha elementos na lista l (lista de queue)
                node, preds, swg = l.pop(0)  # isolar o 1º no na queue
                bwg = 999999
                for elem in self.graph[node]:  # correr os nos no node a pesquisar
                    nwnode, wg = elem
                    if nwnode == d:
                        return preds + [(node, nwnode)] , swg+wg
                    if wg < bwg:
                        bwg = wg
                        nxnode = nwnode
                if nxnode not in vis and nxnode not in l and nxnode != node:  # adicionar á queue
                    l.append((nxnode, [preds + (node, nxnode)], swg + bwg))
                    vis.append(node)
            return None

    def reachable_with_dist(self,s):  # travessia total do grafo mas com as distancias associadas, faz a travessia sobre todos os pontos
        res = []
        l = [(s, 0)]
        while len(l) > 0:
            node, swg = l.pop(0)
            if node != s: res.append((node, swg))
            for elem in self.graph[node]:
                nwnode, wg = elem
                if not is_in_tuple_list(l, nwnode) and not is_in_tuple_list(res,nwnode):  # juntar sempre a distancia ao registro dos elementos nos grafos
                    l.append((nwnode, swg + wg))
        return res

    ## cycles
    def node_has_cycle(self, v):
        l = [v]
        res = False
        visited = [v]
        while len(l) > 0:
            node = l.pop(0)
            for elem in self.graph[node]:
                nwnode,wg = elem
                if nwnode == v:
                    return True
                elif nwnode not in visited:
                    l.append(nwnode)
                    visited.append(nwnode)
        return res

    def has_cycle(self):
        res = False
        for v in self.graph.keys():
            if self.node_has_cycle(v): return True
        return res


def is_in_tuple_list(tl, val):
    res = False
    for (x, y) in tl:
        if val == x: return True
    return res


def test1():
    gr = MyGraph({1: [(2,2)], 2: [(3,4)], 3: [(2,3), (4,2)], 4: [(2,5)]})  # criar o grafo
    gr.print_graph()
    print(gr.get_nodes())
    print(gr.get_edges())


def test2():
    gr2 = MyGraph()
    gr2.add_vertex(1)
    gr2.add_vertex(2)
    gr2.add_vertex(3)
    gr2.add_vertex(4)

    gr2.add_edge(1, 2, 2)
    gr2.add_edge(2, 3, 4)
    gr2.add_edge(3, 2, 3)
    gr2.add_edge(3, 4, 2)
    gr2.add_edge(4, 2, 5)

    gr2.print_graph()


def test3():
    gr = MyGraph({1: [(2,2)], 2: [(3,4)], 3: [(2,3), (4,2)], 4: [(2,5)]})
    gr.print_graph()
    print()
    print(gr.get_successors(2))
    print()
    print(gr.get_predecessors(2))
    print()
    print(gr.get_adjacents(2))
    print()
    print(gr.in_degree(2))
    print()
    print(gr.out_degree(2))
    print()
    print(gr.degree(2))


def test4():
    gr = MyGraph({1: [(2,2)], 2: [(3,4)], 3: [(2,3), (4,2)], 4: [(2,5)]})

    print(gr.distance(1, 4))
    print(gr.distance(4, 3))

    print(gr.shortest_path(1, 4))
    print(gr.shortest_path(4, 3))

    print(gr.reachable_with_dist(1))
    print(gr.reachable_with_dist(3))

    gr2 = MyGraph({1: [(2,2)], 2: [(3,4)], 3: [(2,3), (4,2)], 4: [(2,5)]})

    print(gr2.distance(2, 1))
    print(gr2.distance(1, 5))

    print(gr2.shortest_path(1, 5))
    print(gr2.shortest_path(2, 1))

    print(gr2.reachable_with_dist(1))
    print(gr2.reachable_with_dist(4))


def test5():
    gr = MyGraph({1: [(2,2)], 2: [(3,4)], 3: [(2,3), (4,2)], 4: [(2,5)]})
    print(gr.node_has_cycle(2))
    print(gr.node_has_cycle(1))
    print(gr.has_cycle())

    gr2 = MyGraph({1: [(2,2)], 2: [(3,4)], 3: [(2,3), (4,2)], 4: [(2,5)]})
    print(gr2.node_has_cycle(1))
    print(gr2.has_cycle())


if __name__ == "__main__":
    #test1()
    #test2()
    #test3()
    #test4()
    test5()
