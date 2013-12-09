# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
from numpy.matlib import zeros
from scipy.linalg import solve
import math

import matplotlib.cm as cm
from matplotlib.pyplot import imshow, show, xticks



# <codecell>

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

# TODO: special handling of small cases (n=2)

# The matrix colums are going to be:
# h0 v0 f1 f2 ... fn v1

# The matrix rows are going to be:
# j0.x j0.y j1.x j1y .. jn.x jn.y
# Define functions for getting the matrix indexes
def eq(joint, axis):
    if axis == 'x':
        return 2 * joint
    elif axis == 'y':
        return 2 * joint + 1

def f(force):
    return force + 1


def build_matrix(n, span, h, C):
    h0_index = 0
    v0_index = 1
    v1_index = 4 * n - 1
    c_index = 4 * n  # amount of equations + 1
    m = zeros((4 * n, 4 * n + 1))

    link_x = 1.0 * span / n
    link_y = 1.0 * h
    link_diag = (link_x ** 2 + link_y ** 2) ** 0.5
    y = link_y / link_diag
    x = link_x / link_diag

    from ipdb import set_trace; set_trace()


    ### Quantities
    j_count = 2 * n       # Number of joints
    j_max = j_count - 1   # Last joint index
    l_count = 4 * n - 3   # Number of links
    l_max = l_count    # Last link index

    if n == 2:
        # Fill j0.x equation
        # h0 = f1 + f2.x
        m[eq(0, 'x'), h0_index] = 1
        m[eq(0, 'x'), f(1)] = -1
        m[eq(0, 'x'), f(2)] = -x

        # Fill j0.y equation
        # v0 = f2.y
        m[eq(0, 'y'), v0_index] = 1
        m[eq(0, 'y'), f(2)] = -y

        # Fill j1.x equation
        # f1 = f4
        m[eq(1, 'x'), f(1)] = 1
        m[eq(1, 'x'), f(4)] = -1

        # Fill j1.y equation
        # -c1 = f3
        m[eq(1, 'y'), f(3)] = 1
        m[eq(1, 'y'), c_index] = C[0]

        # Fill j2.x equation
        # f2.x = f5.x + f6
        m[eq(2, 'x'), f(2)] = x
        m[eq(2, 'x'), f(5)] = -x

        # Fill j2.y equation
        # f2.y + f3 + f5.y = 0
        m[eq(2, 'y'), f(2)] = y
        m[eq(2, 'y'), f(3)] = 1
        m[eq(2, 'y'), f(5)] = y

        # Fill j(j_max).x equation
        # f(l_max).x + f(l_max - 1) = 0
        m[eq(j_max, 'x'), f(l_max)] = x
        m[eq(j_max, 'x'), f(l_max - 1)] = 1

        # Fill j(j_max).y equation
        # v1 = f(l_max).y
        m[eq(j_max, 'y'), v1_index] = 1
        m[eq(j_max, 'y'), f(l_max)] = -y

        return m
        pass
    # Fill j0.x equation
    # h0 = f1 + f2.x
    m[eq(0, 'x'), h0_index] = 1
    m[eq(0, 'x'), f(1)] = -1
    m[eq(0, 'x'), f(2)] = -x

    # Fill j0.y equation
    # v0 = f2.y
    m[eq(0, 'y'), v0_index] = 1
    m[eq(0, 'y'), f(2)] = -y

    # Fill j1.x equation
    # f1 = f4
    m[eq(1, 'x'), f(1)] = 1
    m[eq(1, 'x'), f(4)] = -1

    # Fill j1.y equation
    # -c1 = f3
    m[eq(1, 'y'), f(3)] = 1
    m[eq(1, 'y'), c_index] = C[0]

    # Fill j2.x equation
    # f2.x = f5.x + f6
    m[eq(2, 'x'), f(2)] = x
    m[eq(2, 'x'), f(5)] = -x
    m[eq(2, 'x'), f(6)] = -1

    # Fill j2.y equation
    # f2.y + f3 + f5.y = 0
    m[eq(2, 'y'), f(2)] = y
    m[eq(2, 'y'), f(3)] = 1
    m[eq(2, 'y'), f(5)] = y

    # Fill j(n - 1).x: center lower.x equation
    # f(2 * n - 4) + f(2 * n - 3).x = f(2 * n) + f(2 * n + 1)
    m[eq(n - 1, 'x'), f(2 * n - 4)] = 1
    m[eq(n - 1, 'x'), f(2 * n - 3)] = x
    m[eq(n - 1, 'x'), f(2 * n)] = -1
    m[eq(n - 1, 'x'), f(2 * n + 1)] = -x

    # Fill j(n - 1).y: center lower.y equation
    # -c(n/2) = f(2 * n - 3).y + f(2 * n - 1) +  f(2 * n + 1).y
    m[eq(n - 1, 'y'), c_index] = C[n / 2 - 1]
    m[eq(n - 1, 'y'), f(2 * n - 3)] = y
    m[eq(n - 1, 'y'), f(2 * n - 1)] = 1
    m[eq(n - 1, 'y'), f(2 * n + 1)] = y

    # Fill jn.x: center upper.x equation
    # f(2 * n - 2) = f(2 * n + 2)
    m[eq(n, 'x'), f(2 * n - 2)] = 1
    m[eq(n, 'x'), f(2 * n + 2)] = -1

    # Fill jn.y: center upper.y equation
    # f(2 * n - 1) = 0
    m[eq(n, 'y'), f(2 * n - 1)] = 1

    # Fill j(j_max-2).x equation
    # f(l_max -5) = f(l_max - 1)
    m[eq(j_max - 2, 'x'), f(l_max - 5)] = 1
    m[eq(j_max - 2, 'x'), f(l_max - 1)] = -1

    # Fill j(j_max-2).y equation
    # -c(n-2) = f(l_max - 2)
    m[eq(j_max - 2, 'y'), f(l_max - 2)] = 1
    m[eq(j_max - 2, 'y'), c_index] = C[n - 2]

    # Fill j(j_max-1).x equation
    # f(l_max -3)+ f(l_max - 4).x = f(l_max).x
    m[eq(j_max - 1, 'x'), f(l_max - 3)] = 1
    m[eq(j_max - 1, 'x'), f(l_max - 4)] = x
    m[eq(j_max - 1, 'x'), f(l_max)] = -x

    # Fill j(j_max-1).y equation
    # f(l_max - 4).y + f(l_max - 2) + f(l_max).y = 0
    m[eq(j_max - 1, 'y'), f(l_max - 4)] = y
    m[eq(j_max - 1, 'y'), f(l_max - 2)] = 1
    m[eq(j_max - 1, 'y'), f(l_max)] = y

    # Fill j(j_max).x equation
    # f(l_max).x + f(l_max - 1) = 0
    m[eq(j_max, 'x'), f(l_max)] = x
    m[eq(j_max, 'x'), f(l_max - 1)] = 1

    # Fill j(j_max).y equation
    # v1 = f(l_max).y
    m[eq(j_max, 'y'), v1_index] = 1
    m[eq(j_max, 'y'), f(l_max)] = -y

    # Fill j(2k) - 1, 1 < k < n / 2
    for k in xrange(3, n - 1, 2):
        # f(2k -2) + f(2k -1).x = f(2k + 2)
        m[eq(k, 'x'), f(2 * k - 2)] = 1
        m[eq(k, 'x'), f(2 * k - 1)] = x
        m[eq(k, 'x'), f(2 * k + 2)] = -1

        # -c(k/2) = f(2k - 1).y + f(2k + 1)
        m[eq(k, 'y'), c_index] = C[k / 2]
        m[eq(k, 'y'), f(2 * k - 1)] = y
        m[eq(k, 'y'), f(2 * k + 1)] = 1

    # Fill j(2k)    , 1 < k < n / 2
    for k in xrange(4, n - 1, 2):
        # f(2k - 2) = f(2k + 2) + f(2k + 1).x
        m[eq(k, 'x'), f(2 * k - 2)] = -1
        m[eq(k, 'x'), f(2 * k + 2)] = 1
        m[eq(k, 'x'), f(2 * k + 1)] = x

        # f(2k - 1) + f(2k + 1).y = 0
        m[eq(k, 'y'), f(2 * k - 1)] = 1
        m[eq(k, 'y'), f(2 * k + 1)] = y

    # Fill j(2k) - 1, n / 2 + 1 < k < n - 1
    for k in range(n + 1, 2 * n - 3, 2):
        # f(2k - 2) = f(2k + 2) + f(2k + 3).x
        m[eq(k, 'x'), f(2 * k - 2)] = -1
        m[eq(k, 'x'), f(2 * k + 2)] = 1
        m[eq(k, 'x'), f(2 * k + 3)] = x

        # -c[k/2] = f(k * 2 + 1) + f(k * 2 + 3).y
        m[eq(k, 'y'), c_index] = C[k / 2]
        m[eq(k, 'y'), f(2 * k + 1)] = 1
        m[eq(k, 'y'), f(2 * k + 3)] = y

    # Fill j(2k)    , n / 2 + 1 < k < n - 1
    for k in range(n + 2, 2 * n - 3, 2):
        # f(2k -3).x + f(2k - 2) = f(2k + 2)
        m[eq(k, 'x'), f(2 * k - 3)] = x
        m[eq(k, 'x'), f(2 * k - 2)] = 1
        m[eq(k, 'x'), f(2 * k + 2)] = -1

        # f(2k -3).y + f(2k -1) = 0
        m[eq(k, 'y'), f(2 * k - 3)] = y
        m[eq(k, 'y'), f(2 * k - 1)] = 1

    return m

# <rawcell>

# <markdowncell>

# Setting the right-hand side of the equations
# --------------------------------------------
#
# Given the list of weights applied at the junctions, we must that into the right-hand
# side of the equations for the matrix

# <codecell>



def solve_problem(A):
    # Use the scipy.linalg solver
    # Make matrix square and use the weights solutions as the B matrix
    # Solve and show the result
    n, m = np.shape(A)
    square = A[:, 0:n]  # Take the last row out
    B = A[:, n]  # Last row

    solution = solve(square, B)
    return solution


def backward_sub(A):
    n, m = np.shape(A)

    solutions = [0] * (n - 1)
    results = A[:, -1]  # The last column of the augmented matrix

    for i in xrange(n, -1, -1):
        solutions[i] = results[i] - (A[i, i:n] * results[i:n]) / A[i, i]

    return solutions


def show_matrix(m):
    from matplotlib.pylab import matshow
    #im = imshow(m, cmap=cm.Reds_r, interpolation='none')
    #label = ['h0', 'v0', 'f1', 'f2', 'f3', 'f4', 'f1', 'f2', 'f3', 'f4', 'f1', 'f2', 'f3', 'f4', 'f1', 'f2', 'f3', 'f4', 'f1', 'f1']
    #show(im)
    matshow(m, fignum=1, cmap=cm.Reds_r, interpolation='none')
    show()


def parse_input():
    import sys
    filename = sys.argv[1]
    data = open(filename, 'r')
    span = float(data.readline())
    h = float(data.readline())
    n = int(float(data.readline()))

    C = [0] * (n - 1)
    for item in xrange(n - 1):
        C[item] = (float(data.readline()))

    return build_matrix(n, span, h, C)


def check_final_force_results(force_i, n, span, h, param_to_increment='span'):
    force_i_results = []
    if param_to_increment == 'span':
        C = [4] * (n - 1)
        for i in xrange(10,span):
            force_i_results.append(solve_problem(build_matrix(n, i, h, C))[f(force_i)])

    elif param_to_increment == 'c_uniform':
        C = [1] * (n - 1)
        for i in xrange(1,100):
            force_i_results.append(solve_problem(build_matrix(n, span, h, C))[f(force_i)])
            C = [ x + 1 for x in C ]

    elif param_to_increment == 'c_rotative_weigth':
        for i in xrange(0, n - 1):
            C = [4] * (n - 1)
            C[i] = 40
            force_i_results.append(solve_problem(build_matrix(n, span, h, C))[f(force_i)])

    print "Results for the force %s incrementing %s: %s" % (force_i, param_to_increment,[ x[0] for x in force_i_results ])


def get_max_stress(m):
    pass

def compute_cost(m, sections, h, span, max_stress):
    x = span/h
    diagonal = math.sqrt(h ** 2 + x ** 2)
    return (diagonal * sections + x * (2 * sections - 2) + h * (sections - 1)) * max_stress


if __name__ == '__main__':
    m = parse_input()
    #m = build_matrix(2, 4, 2, [1])
    #show_matrix(m)
    dim_n = m.shape[0]
    square = m[:, 0:dim_n]  # Take the last row out
    B = m[:, dim_n]  # Last row
    # from ipdb import set_trace; set_trace()
    #show_matrix(solution)
    #dim_n = np.shape(m)[0]
    # show_matrix(m)
    # check_matrix(m)
    #experimento = BaseExperimento([[18, 2, 6, 1, 1, 1, 1, 1]])
    #NHistStudy()
    result = solve_problem(m)
    #cost =
    #total_cost = cost * (4 + dim_n - 3)


    from experiments import check_matrix, format_result
    check_matrix(m)
    #print total_cost
    print "Results:"
    format_result(result)

# check_final_force_results(9, 6, 100, 4)

# check_final_force_results(9, 6, 100, 4, 'c_uniform')

#check_final_force_results(9, 6, 100, 4, 'c_rotative_weigth')
