# -*- coding: utf-8 -*-

class SuffixTree:
    
    def __init__(self):
        self.nodes = { 0:(-1,{}) } # root node
        self.num = 0
    
    def print_tree(self):
        for k in self.nodes.keys():
            if self.nodes[k][0] < 0:
                print (k, "->", self.nodes[k][1]) 
            else:
                print (k, ":", self.nodes[k][0])
                
    def add_node(self, origin, symbol, leafnum = -1):
        self.num += 1
        self.nodes[origin][1][symbol] = self.num
        self.nodes[self.num] = (leafnum,{})
        
    def add_suffix(self, p, sufnum):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node][1].keys():
                if pos == len(p)-1:
                    self.add_node(node, p[pos], sufnum)
                else:
                    self.add_node(node, p[pos])     
            node = self.nodes[node][1][p[pos]]
            pos += 1
    
    def suffix_tree_from_seq(self, text):
        t = text+"$"
        for i in range(len(t)):
            self.add_suffix(t[i:], i)
            
    def find_pattern(self, pattern):
        pos = 0
        node = 0
        for pos in range(len(pattern)):
            if pattern[pos] in self.nodes[node][1].keys():
                node = self.nodes[node][1][pattern[pos]]
                pos += 1
            else: return None
        return self.get_leafes_below(node)
        

    def get_leafes_below(self, node):
        res = []
        if self.nodes[node][0] >=0: 
            res.append(self.nodes[node][0])            
        else:
            for k in self.nodes[node][1].keys():
                newnode = self.nodes[node][1][k]
                leafes = self.get_leafes_below(newnode)
                res.extend(leafes)
        return res

    def get_nodes_below(self, node):
        res = []
        if node in self.nodes[node][0]:
            res.append(node)
            for n in res:
                if n in self.nodes[n][0]:
                    for k in self.nodes[n][1].keys():
                        newnode = self.nodes[n][1][k]
                        res.append(newnode)
        return res






def test1():
    seq = "TACTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    st.print_tree()
    print(st.get_nodes_below(0))

#def test():
    seq = "TACTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    st.print_tree()
    print (st.find_pattern("TA"))
    print (st.find_pattern("ACG"))

#def test2():
    seq = "TACTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    print (st.find_pattern("TA"))
    #print(st.repeats(2,2))

test1()
#test()
#print()
#test2()
        
            
    
    
