# -*- coding: utf-8 -*-

from Bio import SeqIO
from graphviz import Digraph
from collections import defaultdict

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1


class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.length = 2
        self.coverage = 0
    
    def calc_coverage(self,c1,c2):
        self.coverage = (c1+c2)/2
        
    def increase_length(self):
        self.length+=1


class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_indx in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            if kmer in self.vertices[next_kmer].in_edges.keys():
                self.vertices[next_kmer].in_edges[kmer][0].increase_length()
                self.vertices[kmer].out_edges[next_kmer][0].increase_length()
            else:
                self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
                self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
        
    def draw_graph(self, full=True):
        
        de_bruin_graph=Digraph(comment='Sequence_graph')
        
        def add_vertices(vert):
            if full==True:
                de_bruin_graph.node(vert.seq)
            elif full==False:
                de_bruin_graph.node(vert.seq, str(vert.coverage))

        
        for v in self.vertices.keys():
            add_vertices(self.vertices[v])
            
            
        for v in self.vertices.keys():
            for e in self.vertices[v].out_edges.keys():

                if full==True:
                    edge_label=self.vertices[v].out_edges[e][0].seq
                    
                elif full==False:
                    edge_label=str(self.vertices[v].out_edges[e][0].length)+', '+str(self.vertices[v].out_edges[e][0].coverage)
                    
                de_bruin_graph.edge(v, e, label=edge_label)
        print(de_bruin_graph)
        self.graph_code=de_bruin_graph
        de_bruin_graph.view()
        
    def write_to_file(self, path):
        with open(path, 'w') as dot_file:
            dot_file.write(self.graph_code.source)



if __name__ == '__main__':
    
    dataset = "path_to_input"

    k = 55
    
    my_graph = Graph(k)
    
    with open(dataset, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            read = str(record.seq)
            my_graph.add_read(read)
#            my_graph.add_read( str(record.reverse_complement().seq) )

    my_graph.calc_init_edge_coverage()
'''    
    for v in my_graph.vertices:
        print('Vertex: {}, coverage: {}'.format(v,my_graph.vertices[v].coverage))
        for e in my_graph.vertices[v].out_edges:
            print('-> Out edge: {}'.format(e))
        for e in my_graph.vertices[v].in_edges:
            print('-> In edge: {}'.format(e))
'''
            
my_graph.draw_graph(full=1)

out_file="path_to_output_to_write_dot_to"

my_graph.write_to_file(out_file)