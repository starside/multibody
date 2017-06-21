from subprocess import Popen, PIPE
import os
import re
from graph_tool.all import *
import numpy as np
from scipy.linalg import toeplitz

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

t = buildTorus(50)
print t
graph_draw(t, vertex_font_size=20, output_size=(600, 600) )