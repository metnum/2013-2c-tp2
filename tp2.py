# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from numpy.matlib import zeros
import matplotlib.cm as cm
from matplotlib.pyplot import imshow, show

# <codecell>

### Input
n = 80                # Number of bridge segments
span = 18            # Length of the bridge
h = 2                # Height
C = [x / 80 for x in range(1, 80)]  # n - 1 weights in the lower joints

### Quantities
j_count = 2 * n       # Number of joints
j_max = j_count - 1   # Last joint index
l_count = 4 * n - 3   # Number of links
l_max = l_count -1    # Last link index

# <markdowncell>

# Force layout
# ------------
# 
# Compute the sines and cosines of theta, depending on the quadrant (denoted by ur, ul, ll, lr).

# <rawcell>

# theta_ur = pi/4
# theta_ul = (3/4) * pi
# theta_ll = (5/4) * pi
# theta_lr = (7/4) * pi
# 
# cos_ur = cos(theta_ur)
# cos_ul = cos(theta_ul)
# cos_ll = cos(theta_ll)
# cos_lr = cos(theta_lr)
# 
# sin_ur = sin(theta_ur)
# sin_ul = sin(theta_ul)
# sin_ll = sin(theta_ll)
# sin_lr = sin(theta_lr)
# 
# cos_ur, cos_ul, cos_ll, cos_lr, sin_ur, sin_ul, sin_ll, sin_lr, 

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
# * *c % 2 == 1*: connections to *f(2c - 2)* (horizontal), *f(2c + 1)* (vertical), *f(2c + 2) (horizontal)*
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

# Create the matrix
# We have 4n equations, and an extra column for the C weights
m = zeros((4 * n, 4 * n + 1))

# TODO: special handling of small cases (n=2)

# The matrix colums are going to be:
# h0 v0 f1 f2 ... fn v1

h0_index = 0
v0_index = 1
v1_index = 4 * n - 1
c_index = 4 * n  # amount of equations + 1

link_x = 1.0 * span / n
link_y = 1.0 * h
link_diag = (link_x ** 2 + link_y ** 2) ** 0.5
x = link_y / link_diag
y = link_x / link_diag

# The matrix rows are going to be:
# j0.x j0.y j1.x j1y .. jn.x jn.y

# Define functions for getting the matrix indexes
def eq(joint, axis):
    if axis == 'x': return 2 * joint
    elif axis == 'y': return 2 * joint + 1
    
    
def f(force):
    return force + 1


# Fill j0.x equation
# h0 = f1 + f2.x
m[eq(0, 'x'),h0_index] = 1
m[eq(0, 'x'),f(1)] = -1
m[eq(0, 'x'),f(2)] = -x

# Fill j0.y equation
# v0 = f2.y
m[eq(0, 'y'),v0_index] = 1
m[eq(0, 'y'),f(2)] = -y

# Fill j1.x equation
# f1 = f4
m[eq(1, 'x'),f(1)] = 1
m[eq(1, 'x'),f(4)] = 1

# Fill j1.y equation
# -c1 = f3
m[eq(1, 'y'),f(3)] = 1
m[eq(1, 'y'),c_index] = C[0]

# Fill j2.x equation
# f2.x = f5.x + f6
m[eq(2, 'x'),f(2)] = x
m[eq(2, 'x'),f(5)] = -x
m[eq(2, 'x'),f(6)] = -1

# Fill j2.y equation
# f2.y + f3 + f5.y = 0
m[eq(2, 'y'),f(2)] = y
m[eq(2, 'y'),f(3)] = 1
m[eq(2, 'y'),f(5)] = y

# Fill j(n - 1).x: center lower.x equation
# f(2 * n - 4) + f(2 * n - 3).x = f(2 * n) + f(2 * n + 1)
m[eq(n - 1, 'x'),f(2 * n - 4)] = 1
m[eq(n - 1, 'x'),f(2 * n - 3)] = x
m[eq(n - 1, 'x'),f(2 * n)] = -1
m[eq(n - 1, 'x'),f(2 * n + 1)] = -x

# Fill j(n - 1).y: center lower.y equation
# -c(n/2) = f(2 * n - 3).y + f(2 * n - 1) +  f(2 * n + 1).y
m[eq(n - 1, 'y'),c_index] = C[n/2 - 1]
m[eq(n - 1, 'y'),f(2 * n - 3)] = y
m[eq(n - 1, 'y'),f(2 * n - 1)] = 1
m[eq(n - 1, 'y'),f(2 * n + 1)] = y

# Fill jn.x: center upper.x equation
# f(2 * n - 2) = f(2 * n + 2)
m[eq(n, 'x'),f(2 * n - 2)] = 1
m[eq(n, 'x'),f(2 * n + 2)] = -1

# Fill jn.y: center upper.y equation
# f(2 * n - 1) = 0
m[eq(n, 'x'),f(2 * n - 1)] = 1

# Fill j(j_max-2).x equation
# f(l_max -5) = f(l_max - 1)
m[eq(j_max - 2, 'x'),f(l_max - 5)] = 1
m[eq(j_max - 2, 'x'),f(l_max - 1)] = -1

# Fill j(j_max-2).y equation
# -c(n-2) = f(l_max - 2)
m[eq(j_max - 2, 'y'),f(l_max - 1)] = 1
m[eq(j_max - 2, 'y'),c_index] = C[n - 2]

# Fill j(j_max-1).x equation
# f(l_max -3)+ f(l_max - 4).x = f(l_max).x
m[eq(j_max - 1, 'x'),f(l_max - 3)] = 1
m[eq(j_max - 1, 'x'),f(l_max - 4)] = x
m[eq(j_max - 1, 'x'),f(l_max)] = -x

# Fill j(j_max-1).y equation
# f(l_max - 4).y + f(l_max - 2) + f(l_max).y = 0
m[eq(j_max - 1, 'y'),f(l_max - 4)] = y
m[eq(j_max - 1, 'y'),f(l_max - 2)] = 1
m[eq(j_max - 1, 'y'),f(l_max)] = y

# Fill j(j_max).x equation
# f(l_max).x + f(l_max - 1) = 0
m[eq(j_max, 'x'),f(l_max)] = x
m[eq(j_max, 'x'),f(l_max - 1)] = 1

# Fill j(j_max).y equation
# v1 = f(l_max).y
m[eq(j_max, 'y'),v1_index] = 1
m[eq(j_max, 'y'),f(l_max)] = -y

# Fill j(2k) - 1, 1 < k < n / 2
for k in xrange(3, n - 1, 2):
    # f(2k -2) + f(2k -1).x = f(2k + 2)
    m[eq(k, 'x'),f(2 * k - 2)] = 1
    m[eq(k, 'x'),f(2 * k - 1)] = x
    m[eq(k, 'x'),f(2 * k + 2)] = -1
    
    # -c(k/2) = f(2k - 1).y + f(2k + 1)
    m[eq(k, 'y'),c_index] = C[k/2]
    m[eq(k, 'y'),f(2 * k - 1)] = y
    m[eq(k, 'y'),f(2 * k + 1)] = 1
    
# Fill j(2k)    , 1 < k < n / 2
for k in xrange(4, n - 1, 2):
    # f(2k - 2) = f(2k + 2) + f(2k + 1).x
    m[eq(k, 'x'),f(2 * k - 2)] = -1
    m[eq(k, 'x'),f(2 * k + 2)] = 1
    m[eq(k, 'x'),f(2 * k + 1)] = x
    
    # f(2k - 1) + f(2k + 1).y = 0
    m[eq(k, 'y'),f(2 * k - 1)] = 1
    m[eq(k, 'y'),f(2 * k + 1)] = y

# Fill j(2k) - 1, n / 2 + 1 < k < n - 1
for k in range(n + 1, 2 * n - 3, 2):
    # f(2k - 2) = f(2k + 2) + f(2k + 3).x
    m[eq(k, 'x'),f(2 * k - 2)] = -1
    m[eq(k, 'x'),f(2 * k + 2)] = 1
    m[eq(k, 'x'),f(2 * k + 3)] = x
    
    # -c[k/2] = f(k * 2 + 1) + f(k * 2 + 3).y
    m[eq(k, 'y'),c_index] = C[k/2]
    m[eq(k, 'y'),f(2 * k + 1)] = 1
    m[eq(k, 'y'),f(2 * k + 3)] = y
                  
# Fill j(2k)    , n / 2 + 1 < k < n - 1
for k in range(n + 2, 2 * n - 3, 2):
    # f(2k -3).x + f(2k - 2) = f(2k + 2)
    m[eq(k, 'x'),f(2 * k - 3)] = 1
    m[eq(k, 'x'),f(2 * k - 2)] = 1
    m[eq(k, 'x'),f(2 * k + 2)] = -1
    
    # f(2k -3).y + f(2k -1) = 0
    m[eq(k, 'y'),f(2 * k - 3)] = y
    m[eq(k, 'y'),f(2 * k - 1)] = 1

# <rawcell>

# # Forces for v0
# m[0, 0] = 0  # No horizontal force
# m[0, 2] = 1  # F1
# m[0, 3] = cos_ur  # Cos(theta) * F2
# m[1, 1] = 1  # Weight force v0
# m[1, 3] = sin_ur  # Sin(theta) * F2
# 
# # Forces for v1
# m[(4*n)-2, 21] = 1
# m[(4*n)-2, 22] = cos_ul
# m[(4*n)-1, 22] = sin_ul
# m[(4*n)-1, 23] = 1

# <rawcell>

# # Set the rest of the matrix
# for c in range(1, 10):
#     for r in range(2*c, 2*c + 2):
#         # Horizontal equation
#         if r % 2 == 0:
#             # Top vertex
#             if c % 2 == 0:
#                 m[r, 2*c - 1] = 1
#                 m[r, 2*c + 3] = 1
#                 if 1 < c < n:  # Diagonal links
#                     m[r, 2*c + 2] = cos_lr
#                 elif n < c < 2*n:
#                     m[r, 2*c - 2] = cos_ll
#             # Bottom vertex
#             else:
#                 m[r, 2*c + 3] = 1
#                 m[r, 2*c - 1] = 1
#                 # First vertex has no upper right
#                 # diagonal link
#                 if 1 < c < n: # Diagonal links
#                     m[r, 2*c] = cos_ul
#                 # Last rightmost vertex does not have any
#                 # block to the right to connect to.
#                 elif n < c < 2*n:
#                     m[r, 2*c + 4] = cos_ur
#         # Vertical equation
#         else:
#             if c % 2 == 0:
#                 m[r, 2*c] = 1
#                 if 1 < c < n:
#                     m[r, 2*c + 2] = sin_lr
#                 elif n < c < 2*n:
#                     m[r, 2*c - 2] = sin_ll
#             else:
#                 m[r, 2*c + 2] = 1
#                 if 1 < c < n:
#                     m[r, 2*c] = sin_ul
#                 elif n < c < 2*n:
#                     m[r, 2*c + 4] = sin_ur

# <markdowncell>

# Setting the right-hand side of the equations
# --------------------------------------------
# 
# Given the list of weights applied at the junctions, we must that into the right-hand
# side of the equations for the matrix

# <codecell>

print m
im = imshow(m, cmap = cm.Reds_r, interpolation='none')
show(im)

# <codecell>


# <codecell>


