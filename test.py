from subprocess import Popen, PIPE
import os
import re
from graph_tool.all import *
import numpy as np
from scipy.linalg import toeplitz

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
    assert(data is not None)
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


options = {"dendrimer expansion": dendrimer_expansion,
           "small example": small_example,
           "linear expansion": linear_expansion,
           "torus expansion":torus_expansion
           }

options["linear expansion"]()
