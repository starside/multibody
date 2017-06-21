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
    for row in range(width):
        for col in range(width):
            m = row*width + col
            if col != width - 1:
                g.add_edge(m, m+1)
            if row != width - 1:
                g.add_edge(m, m + width)
    
    return g

t = buildTorus(4)
print t