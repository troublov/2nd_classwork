# -*- coding: utf-8 -*-
"""
@author: orlov
"""

from Bio import SeqIO
from graphviz import Digraph
from collections import defaultdict

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        self.been_here=False
        
    def increase_coverage(self):
        self.coverage += 1
    
    def merge(self, next_one):
        self.seq=self.seq+next_one[-1]
        
    def was_here(self):
        self.been_here=True

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

        self.graph_code=de_bruin_graph
        de_bruin_graph.view()
        
    def write_to_file(self, path):
        with open(path, 'w') as dot_file:
            dot_file.write(self.graph_code.source)
            
        

            
    def compress_graph(self):
        
        def contig_start_search(vrt, start, first_vrt):
            its_start_point=start
            
            if len(self.vertices[vrt].in_edges)==1:
                prev_vrt=list(self.vertices[vrt].in_edges)[0]
                
                if self.vertices[prev_vrt].been_here==False:
                    if len(self.vertices[prev_vrt].out_edges)==1 and len(self.vertices[prev_vrt].in_edges)==1:
                        self.vertices[vrt].been_here=True
                        contig_start_search(prev_vrt, False, first_vrt)
                    
                    elif ((len(self.vertices[prev_vrt].out_edges)>1) or (len(self.vertices[prev_vrt].in_edges)>1)) and its_start_point==False and (len(self.vertices[list(self.vertices[vrt].out_edges)[0]].out_edges)==1):
                        self.vertices[vrt].been_here=True
                        starting_vertices_to_compress.append(vrt)
                
                elif ((len(self.vertices[prev_vrt].out_edges)>1) or (len(self.vertices[prev_vrt].in_edges)>1)) and its_start_point==False:
                    self.vertices[vrt].been_here=True
                    starting_vertices_to_compress.append(vrt)
                    
                elif first_vrt==prev_vrt:##if we started from the begining of cycle
                    self.vertices[vrt].been_here=True
                    starting_vertices_to_compress.append(vrt)
                
                
            elif len(self.vertices[vrt].in_edges)==0 and len(self.vertices[vrt].out_edges)==1:
                next_vrt=list(self.vertices[vrt].out_edges)[0]
                
                if len(self.vertices[next_vrt].in_edges)==1 and len(self.vertices[next_vrt].out_edges)==1:
                    self.vertices[vrt].been_here=True
                    starting_vertices_to_compress.append(vrt)
        
                    
        def compress_contig(vrt):
            current_vrt=vrt
            continue_compression=True
            
            while continue_compression==True:
                next_vrt=list(self.vertices[current_vrt].out_edges)[0]
                
                if len(self.vertices[next_vrt].in_edges)>1:
                    continue_compression=False
                    break
                
                elif (len(list(self.vertices[current_vrt].in_edges))!=0 and (next_vrt==list(self.vertices[current_vrt].in_edges)[0])):
                    continue_compression=False
                    break
                
                new_vrt_seq=current_vrt+next_vrt[-1]
                new_vrt_coverage=self.vertices[current_vrt].coverage+self.vertices[next_vrt].coverage
                
                if new_vrt_seq in self.vertices.keys():
                    self.vertices[new_vrt_seq].coverage+=new_vrt_coverage
                    
                else:
                    self.vertices[new_vrt_seq]=Vertex(new_vrt_seq)
                    self.vertices[new_vrt_seq].coverage=new_vrt_coverage
                
                for out_edge_vertex in self.vertices[next_vrt].out_edges.keys():
                    new_out_edge=Edge(new_vrt_seq, out_edge_vertex)
                    
                    self.vertices[new_vrt_seq].out_edges[out_edge_vertex]=[new_out_edge]
                    self.vertices[out_edge_vertex].in_edges[new_vrt_seq]=[new_out_edge]
                    del self.vertices[out_edge_vertex].in_edges[next_vrt]
                    
                for in_edge_vertex in self.vertices[current_vrt].in_edges.keys():
                    new_in_edge=Edge(in_edge_vertex, new_vrt_seq)
                    new_in_edge.seq=in_edge_vertex[0]+new_vrt_seq
                    
                    self.vertices[new_vrt_seq].in_edges[in_edge_vertex]=[new_in_edge]
                    self.vertices[in_edge_vertex].out_edges[new_vrt_seq]=[new_in_edge]
                    del self.vertices[in_edge_vertex].out_edges[current_vrt]
                    
                del self.vertices[current_vrt]
                del self.vertices[next_vrt]
                    
                current_vrt=new_vrt_seq
                
                if len(self.vertices[current_vrt].out_edges)>1 or len(self.vertices[current_vrt].out_edges)==0:
                    continue_compression=False
                
        starting_vertices_to_compress=[]
        
        for vertex in self.vertices.keys():
            if self.vertices[vertex].been_here==False:
                contig_start_search(vertex, True, vertex)
                
                
        for vertex in starting_vertices_to_compress:
            compress_contig(vertex)




if __name__ == '__main__':
    
    dataset = "path_to_input"

    k = 10
    
    my_graph = Graph(k)
    
    with open(dataset, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            read = str(record.seq)
            my_graph.add_read(read)

    my_graph.calc_init_edge_coverage()


    my_graph.draw_graph(full=1)
    my_graph.compress_graph()        
    my_graph.draw_graph(full=1)
    
    
    my_graph.calc_init_edge_coverage()

