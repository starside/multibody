from subprocess import Popen, PIPE
import os
import re
from graph_tool.all import *
import numpy as np
from scipy.linalg import toeplitz
import random

def buildDendrimer(gen, func):
    g = Graph(directed=False)
    v0 = g.add_vertex()  # 0th generation
    vlist = []  # create a list of the external vertices
    # initialize
    for j in range(func):  # 1st generation
        nv1 = g.add_vertex()
        g.add_edge(v0, nv1)
        vlist.append(nv1)
    i = 1  # number of generations
    while i < gen:  # start adding vertices
        terminals = []
        for rv in vlist:
            for j in range(func - 1):
                # add the new vertex
                nv1 = g.add_vertex()
                g.add_edge(rv, nv1)
                terminals.append(nv1)
        vlist = terminals
        i = i + 1
    return g

def buildTorus(width):
    g = Graph(directed=False)
    for i in range(width*width):    # Add vertices
        g.add_vertex()
    #Build grid mesh
    for row in range(width):
        for col in range(width):
            m = row*width + col
            if col != width - 1:
                g.add_edge(m, m+1)
            if row != width - 1:
                g.add_edge(m, m + width)
    #Connect edges in period boundary condition
    for row in range(width):
        g.add_edge(row*width, row*width + width - 1)
    for col in range(width):
        g.add_edge(col, width*width - width + col)
    return g


def buildLinear(N):
	x = np.zeros((N,1))
	x[1] = 1
	return toeplitz(x)


def calcCoeffs(adj, size, dim, eps, a):
    inpt = str(size) + " " + str(dim) + " " + \
        str(eps) + " " + str(a) + " " + adj + "\n"
    command = os.getcwd() + "/multibody"
    p = Popen(command, stdin=PIPE, stdout=PIPE)
    res = p.communicate(input=inpt)
    # Parse with regular expression
    pattern = """^	#Beginning
				(-*\d*\.*\d*(?:e\+){0,1}\d*)\s+	#First number
				(-*\d*\.*\d*(?:e\+){0,1}\d*)\s+	#Second numbert
				(-*\d*\.*\d*(?:e\+){0,1}\d*)\s*		#Third number"""
    expr = re.compile(pattern, re.VERBOSE)
    data = expr.match(res[0])
    try:
        assert(data is not None)
    except:
        raise AssertionError
    # Convert extracted data to floating point
    [rg2, expansion, thirdorder] = [float(i) for i in data.groups()]
    return [rg2, expansion, thirdorder]


def adjToString(matrix):
    matrix = np.array(matrix).reshape(-1)
    string = ""
    for i in range(len(matrix)):
        string += str(matrix[i]) + " "
    string += "\n"
    return string


def linear_expansion():
    for N in range(20, 2000, 30):
        adj = buildLinear(N)
        adjstring = adjToString(adj)
        confinement_energy = 5.0 #Totan energy of confinement
        D = 3
        eps = 0
        a = 1.0
        coefs = calcCoeffs(adjstring, N, D, eps, a)
        #2*std::pow(2.0*M_PI*a*a / D, -D / 2.0)
        z = 2*np.power(2.0*np.pi*a*a / D, -D / 2.0)*np.power(N, 0.5)
        print("{0} {1} {3} {2}".format(
            N * N * np.power(coefs[0], -D / 2.0), coefs[1]/z, N, coefs[2]/(z*z) ))

def dendrimer_expansion():
    genmax = 10
    func = 3
    for gen in range(1, genmax):
        dg = buildDendrimer(gen, func)
        adj = adjacency(dg).todense()
        adjstring = adjToString(adj)
        N = len(adj)
        D = 3
        eps = 10
        a = 1.0
        coefs = calcCoeffs(adjstring, N, D, eps, a)
        print("{0} {1} {2}".format(
            N * N * np.power(coefs[0], -D / 2.0), coefs[1], N))

def torus_expansion():
    sizemin = 4
    sizemax = 50
    for w in range(sizemin, sizemax, 2):
        dg = buildTorus(w)
        adj = adjacency(dg).todense()
        adjstring = adjToString(adj)
        N = len(adj)
        D = 3
        eps = 0
        a = 1.0
        coefs = calcCoeffs(adjstring, N, D, eps, a)
        print("{0} {1} {2}".format(
            N * N * np.power(coefs[0], -D / 2.0), coefs[1], N))

def count_edges(g):
    #num_edges is wrong
    count = 0
    for e in g.edges():
        count = count + 1
    return count

def removeRandomEdge(g):
    """
    removes a random edge from graph g.  Not efficient in implementation.  If a disconnected
    fragment is created, we delete the smallest fragment.  If no fragment was produced, return
    None
    """
    class Frag(BFSVisitor):
        def __init__(self):
            self.vertexmap = dict()

        def discover_vertex(self, u):
            self.vertexmap[u] = None

        def graph(self, g):
            # Build map from old vertices, 
            # to vertices in new graph (starting at 0)
            for i,v in enumerate(self.vertexmap.iterkeys()):
                self.vertexmap[v] = i
            # Create empty graph
            newg = Graph(directed=False)
            #Add Edges
            if len(self.vertexmap) == 1:
                newg.add_vertex()
            else:
                for e in g.edges():
                    # Add edges only in the fragment.  Not efficient.  I should not use graph-tool
                    if e.source() in self.vertexmap and e.target() in self.vertexmap:
                        newg.add_edge(self.vertexmap[e.source()], self.vertexmap[e.target()] )
            return newg


    def countFragment(g, node):
        """ Counts the size of a fragment"""
        fragment = Frag()
        bfs_search(g, node, fragment)
        fragment.size = len(fragment.vertexmap)
        return fragment
    
    numedges = g.num_edges()
    inf = g.num_vertices()*2

    small_fragment = None
    randomedge = random.randint(0, count_edges(g) - 1)
    for i,e in enumerate(g.edges()): #Not efficient way to do this
        if i == randomedge:
            [s, t] = e.source(), e.target()
            g.remove_edge(e)
            fragment1 = countFragment(g, s)
            fragment2 = countFragment(g, t) 
            # delete smallest fragment
            if fragment1.size != g.num_vertices():
                todie = fragment2
                if fragment1.size < fragment2.size:
                    todie =fragment1
                #print len(fragment1.edges), len(fragment2.edges), len(todie.edges)
                #graph_draw(g, vertex_font_size=10, output_size=(800, 800))
                small_fragment = todie.graph(g)
                # remove left over edges that do not get deleted with vertices
                killedges = []
                for deadedge in g.edges():
                    if todie.vertexmap.has_key(deadedge.source()) or todie.vertexmap.has_key(deadedge.target()):
                        killedges.append(deadedge)
                for deadedge in killedges:
                    g.remove_edge(deadedge)
                # Delete the vertices in the smallest fragment
                g.remove_vertex([ dk for dk in todie.vertexmap.iterkeys() ])
            break
    return small_fragment

def torus_percolation():
    gridsize = 32

    dg = buildTorus(gridsize)
    adj = adjacency(dg).todense()
    adjstring = adjToString(adj)
    N = len(adj)
    fragments = []
    print("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}".format(
        "R_G^2", "alpha^2", "N (big fragment)", "# edges (b.f)", "avg. R_g^2 (s.f)", \
        "avg. alpha^2 (s.f)", "avg. mass (s.f)", "# small fragments" ))
    while N > 4:
        N = len(adj)
        D = 3
        eps = 0.0
        a = 1.0
        # Calculate the Rg2, alpha, etc for the large fragment
        coefs = calcCoeffs(adjstring, N, D, eps, a)
        # Delete a random edge, and add the small fragment to a list
        _ = removeRandomEdge(dg)
        if _ is not None:
            fragments.append(_)

        # Calculate avergage R_g, \alpha and fragment mass
        average_rg2 = 0
        average_alpha = 0
        average_mass = 0
        if len(fragments) > 0:
            for f in fragments:
                if f.num_vertices() > 1:
                    # Calculate coefficients of each small fragment
                    [frag_rg2, frag_alpha, _] = calcCoeffs(adjToString(adjacency(f).todense()), f.num_vertices(), D, eps, a)
                    average_rg2 += frag_rg2
                    average_alpha += frag_alpha
                average_mass += f.num_vertices()
            average_alpha = average_alpha / (1.0* len(fragments))
            average_rg2 = average_rg2 / (1.0* len(fragments))
            average_mass = average_mass / (1.0* len(fragments))
            # Calculate average fragment mass

        # Print the results of each calculation
        print("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}".format(
            coefs[0], coefs[1], N, count_edges(dg), average_rg2, average_alpha, average_mass, len(fragments) ))
        # rebuild adjaceny matrix
        adj = adjacency(dg).todense()
        adjstring = adjToString(adj)


def small_example():
    gas = False
    N = 36
    D = 3
    epsmin = 0.01
    epsmax = 20.0
    runs = 100
    a = 1.0

    adj = r"0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0"

    if gas:
        adj = "0 " * (N * N)

    for i in range(runs):
        coefs = calcCoeffs(
            adj, N, D, ((epsmax - epsmin) / runs) * i + epsmin, a)
        print("{0} {1} {2}".format(coefs[0], coefs[1], coefs[2]))

def test():
    dg = buildDendrimer(2,3)
    dg.add_edge(1,4)
    dg.add_edge(1,5)
    dg.add_edge(2,6)
    dg.add_edge(2,7)
    dg.add_edge(3,8)
    dg.add_edge(3,9)
    for i in range(20):
        print "Num edges ", count_edges(dg)
        for e in dg.edges():
            print e
        ng = removeRandomEdge(dg)
        if ng is not None:
            print "Small chunk:", ng.num_vertices()
            for e in ng.edges():
                print e
        print "*****************************"


options = {"dendrimer expansion": dendrimer_expansion,
           "small example": small_example,
           "linear expansion": linear_expansion,
           "torus expansion":torus_expansion,
           "torus percolation":torus_percolation,
           "test":test
           }

options["torus percolation"]()
