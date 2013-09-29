# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from numpy.matlib import matrix, zeros
from math import sin, cos, pi
import matplotlib.cm as cm
from matplotlib.pyplot import imshow, show

# <codecell>

n = 6  # Number of bridge segments
h = None # Height. TODO: Implement
cs = 2*n - 2

l = 4*n - 3   # Number of links
m = zeros((6*4, 6*4 + 1))

# <codecell>

# Create array of weight values

weights = matrix([0, 0, 0, 2, 0, 0, 0]).transpose()

# <markdowncell>

# Force layout
# ------------
#
# Compute the sines and cosines of theta, depending on the quadrant (denoted by ur, ul, ll, lr).

# <codecell>

theta_ur = pi/4
theta_ul = (3/4) * pi
theta_ll = (5/4) * pi
theta_lr = (7/4) * pi

cos_ur = cos(theta_ur)
cos_ul = cos(theta_ul)
cos_ll = cos(theta_ll)
cos_lr = cos(theta_lr)

sin_ur = sin(theta_ur)
sin_ul = sin(theta_ul)
sin_ll = sin(theta_ll)
sin_lr = sin(theta_lr)

# <markdowncell>

# Special vertexes
# ---------------
#
# Define equations for special vertexes v0 and v1
#
# v0:
# h0 = f1 + cos_ur * f2
#
# v0 = sin_ur * f2
#
# v1:
# 0 = f20 + cos_ul * f21
#
# v1 = sin_ul * f21

# <markdowncell>

# Characterizing links
# --------------------
#
# Given 4n - 3 links, any links before the (2n-1)th link are on the left of the bridge, the (2n-1)th is in the middle, and any afterwards are located to the right.
#
# Links whose number is:
#
# * *n % 4 == 0, and n == 1*: located on the bottom of the bridge.
# * *n % 4 == 1, and n == 2*, are the diagonal bridge. If the link number is less than the middle link number, then they each go from the upper left to the lower right diagonals, and from the lower left to the upper right otherwise.
# * *n % 4 == 3* are vertical links.
# * *n % 4 == 2*, except for n == 2, are the top horizontal links
#
# Relation between links and vertices
# -----------------------------------
#
# For the v0 and v1 vertices the connections have been described in the prior sections..
# For the vertices in the squares, we have the following, where n is the number of sections on the bridge:
#
# * *c % 2 == 1*: joint is on the bottom
# * *c % 2 == 0*: joint is on top
# * *c % 2 == 1: connections to *f(2c - 2)* (horizontal), *f(2c + 1)* (vertical), *f(2c + 2) (horizontal)*
# * *c % 2 == 1 and 1 < c < n*: additional connection to *f(2c - 1)* (horizontal and vertical)
# * *c % 2 == 1 and n < c < 2n*: additional connection to *f(2c + 3)* (horizontal and vertical)
# * *c % 2 == 0*: connections to *f(2c - 2)* (horizontal), *f(2c - 1)* (vertical), *f(2c + 2)* (horizonal)
# * *c % 2 == 0 and 1 < c < n*: additional connection to *f(2c + 1)* (horizonal and vertical)
# * *c % 2 == 0 and n < c < 2n*: additional connection to *f(2c - 3)* (horizontal and vertical)

# <markdowncell>

# Matrix layout
# -------------
#
# The matrix will be laid out as such:
#
# * The first column will provide the multiplier for h0
# * The second column will provide the multiplier for v0
# * The last column will provide the multiplier for v1
# * Each column in between v0 and v1 will show Force(X), where X = n -1
# * Every even row will describe an equation for horizontal forces
# * Every odd row will describe an equation for vertical forces

# <markdowncell>

# Filling initial system
# ----------------------
#
# We must fill the spaces in the original matrix with the corresponding data.
# We do so for the special vertices *v0* and *v1*.

# <codecell>

# Forces for v0
m[0, 0] = 0  # No horizontal force
m[0, 2] = 1  # F1
m[0, 3] = cos_ur  # Cos(theta) * F2
m[1, 1] = 1  # Weight force v0
m[1, 3] = sin_ur  # Sin(theta) * F2

# Forces for v1
m[(4*n)-2, 21] = 1
m[(4*n)-2, 22] = cos_ul
m[(4*n)-1, 22] = sin_ul
m[(4*n)-1, 23] = 1

# <codecell>

# Set the rest of the matrix
for c in range(1, 10):
    for r in range(2*c, 2*c + 2):
        # Horizontal equation
        if r % 2 == 0:
            # Top vertex
            if c % 2 == 0:
                m[r, 2*c - 1] = 1
                m[r, 2*c + 3] = 1
                if 1 < c < n:  # Diagonal links
                    m[r, 2*c + 2] = cos_lr
                elif n < c < 2*n:
                    m[r, 2*c - 2] = cos_ll
            # Bottom vertex
            else:
                m[r, 2*c + 3] = 1
                m[r, 2*c - 1] = 1
                # First vertex has no upper right
                # diagonal link
                if 1 < c < n: # Diagonal links
                    m[r, 2*c] = cos_ul
                # Last rightmost vertex does not have any
                # block to the right to connect to.
                elif n < c < 2*n:
                    m[r, 2*c + 4] = cos_ur
        # Vertical equation
        else:
            if c % 2 == 0:
                m[r, 2*c] = 1
                if 1 < c < n:
                    m[r, 2*c + 2] = sin_lr
                elif n < c < 2*n:
                    m[r, 2*c - 2] = sin_ll
            else:
                m[r, 2*c + 2] = 1
                if 1 < c < n:
                    m[r, 2*c] = sin_ul
                elif n < c < 2*n:
                    m[r, 2*c + 4] = sin_ur


# <markdowncell>

# Setting the right-hand side of the equations
# --------------------------------------------
#
# Given the list of weights applied at the junctions, we must that into the right-hand
# side of the equations for the matrix

# <codecell>
im = imshow(m, cmap = cm.Reds_r, interpolation='none')
show(im)

# <markdowncell>

# Graph Plot of bridge forces
# ---------------------------
#
# We use networkx to plot the bridge diagram
import networkx as nx

plot_matrix = m[::2]
G = nx.from_numpy_matrix(m)
pos = {'v0': (0, 0),
       'v1': (n, 0)}

for vertex in range(1, cs + 1):
    coords = (coord_x, coord_y) = (
        (1 + vertex) / 2,
        (vertex - 1) % 2)

    pos[vertex + 1] = coords

nx.draw_networkx_nodes(


