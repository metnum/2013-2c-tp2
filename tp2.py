# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
from numpy.matlib import zeros
from scipy.linalg import solve

import matplotlib.cm as cm
from matplotlib import pyplot as plt
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
    x = link_y / link_diag
    y = link_x / link_diag


    ### Quantities
    j_count = 2 * n       # Number of joints
    j_max = j_count - 1   # Last joint index
    l_count = 4 * n - 3   # Number of links
    l_max = l_count    # Last link index

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


def lu( A ):
    """Factor A into LU by Gaussian elimination with scaled partial pivoting

    USAGE:
        p = lu( A )

    INPUT:
        A     - NumPy matrix to be factored.  This matrix should should
                have a floating point data type.
                ***** THIS MATRIX IS MODIFIED *****

    OUTPUT:
        p     - (integer list) contains the row permutation information

    NOTES:
        This function performs Gaussian elimination with scaled partial
        pivoting on a square matrix A.  It returns the LU factorization of a
        and a permutation vector that indicates the true position of the rows
        in the result.

        The companion function solve() can be used to solve Ax=b once A has
        been factored by this function.
    """

    n, m = np.shape( A )
    if n != m:
        print "Error: input matrix is not square"
        return None

    # Generate initial index vector

    p = range( n )

    # Determine the largest (in magnitude) element in each row.  These
    # factors are used to scale the pivot elements for comparison purposes
    # when deciding which row to use as a pivot row.

    s = [0] * n
    for i in xrange( n ):
        smax = 0.0
        for j in xrange( n ):
            smax = max( smax, abs( A[i,j] ) )
        s[i] = smax

    # Begin Gaussian elimination.

    for k in xrange( n - 1 ):

        # Find the remaining row with the largest scaled pivot.

        rmax = 0.0
        for i in xrange( k, n ):
            print s[p[i]]
            r = abs( A[p[i],k] / s[p[i]] )
            if r > rmax:
                rmax = r
                j = i

        # Row j has the largest scaled pivot, so "swap" that row with the
        # current row (row k).  The swap is not actually done by copying rows,
        # but by swaping two entries in an index vector.

        p[j], p[k] = ( p[k], p[j] )

        # Now carry out the next elimination step as usual, except for the
        # added complication of the index vector.

        for i in xrange( k + 1, n ):
            xmult = A[p[i],k] / A[p[k],k]
            A[p[i],k] = xmult
            for j in xrange( k + 1, n ):
                A[p[i],j] = A[p[i],j] - xmult * A[p[k],j]

    # All done, return factored matrix A and permutation vector p

    return ( A, p )


def solve_lu( A, p, b ):
    """Solves Ax = b given an LU factored matrix A and permuation vector p

    USAGE:
        x = solve( A, p, b )

    INPUT:
        A     - NumPy matrix that contains LU factored form of A
        p     - list or NumPy array of permuted row positions
        b     - NumPy matrix or array containing RHS vector

    OUTPUT:
        x     - (NumPy array) Solution vector

    NOTES:
        This function performs the backward and forward solves necessary to
        use the LU factorization of a square matrix A to solve Ax=b.  The LU
        factorization is assumed to have been generated by the function
        lu() and has an associated permutation vector that indicates the
        true row values in the factored matrix.
    """

    n, m = np.shape( A )
    if n != m:
        print "Error: input matrix is not square"
        return None

    # Forward solve

    x = np.zeros( n )

    for k in xrange( n - 1 ):
        for i in xrange( k + 1, n ):
            b[p[i]] = b[p[i]] - A[p[i],k] * b[p[k]]

    # Backward solve

    for i in xrange( n - 1, -1, -1 ):
        sum = b[p[i]]
        for j in xrange( i + 1, n ):
            sum = sum - A[p[i],j] * x[j]
        x[i] = sum / A[p[i],i]

    # All done, return solution vector

    return x


def gauss(A, b):
    '''
    Gaussian elimination with partial pivoting.
    % input: A is an n x n nonsingular matrix
    %        b is an n x 1 vector
    % output: x is the solution of Ax=b.
    % post-condition: A and b have been modified.
    '''

    n =  len(A)
    if b.size != n:
        raise ValueError("Invalid argument: incompatible sizes between A & b.", b.size, n)
    # k represents the current pivot row. Since GE traverses the matrix in the upper
    # right triangle, we also use k for indicating the k-th diagonal column index.
    for k in xrange(n-1):
        #Choose largest pivot element below (and including) k
        maxindex = abs(A[k:,k]).argmax() + k
        if A[maxindex, k] == 0:
            raise ValueError("Matrix is singular.")
        #Swap rows
        if maxindex != k:
            A[[k,maxindex]] = A[[maxindex, k]]
            b[[k,maxindex]] = b[[maxindex, k]]
        for row in xrange(k+1, n):
            multiplier = A[row, k]/A[k, k]
            #the only one in this column since the rest are zero
            A[row, k] = multiplier
            for col in xrange(k + 1, n):
                A[row, col] = A[row, col] - multiplier*A[k, col]
            #Equation solution column
            b[row] = b[row] - multiplier*b[k]

    x = np.zeros(n)
    k = n-1
    x[k] = b[k]/A[k,k]
    while k >= 0:
        x[k] = (b[k] - np.dot(A[k,k+1:],x[k+1:]))/A[k,k]
        k = k-1
    return x


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


def force_position(pos):
    if pos == 0:
        return "H0"
    if pos == 1:
        return "V0"
    if pos < 4 * n - 1:
        return "F" + str(pos - 1)
    else:
        return "V1"


def format_result(result):
    """
    Given a result column vector, show the values assigned
    """
    for pos in range(0, result.shape[0]):
        print force_position(pos) + ": " + str(result[pos, 0])


def check_matrix(m):
    # Horizontal and vertical link columns should add to 2
    # Diagonal links should add to 2/sqrt(2)

    def print_res(r):
        for result in r:
            print result[0], result[1]

    print "Getting sums of lower horizonal links"
    n = np.shape(m)[0] / 4
    sums = []
    for force in range(4, 4*n - 3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:,f(force)]]))
    print_res(zip(range(4, 4*n, 4), sums))


    print "Getting sums of upper horizonal links"
    sums = []
    for force in range(6, 4*n - 3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:,f(force)]]))
    print_res(zip(range(6, 4*n, 4), sums))


    print "Getting sums of vertical links"
    sums = []
    for force in range(3, 4*n -3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:,f(force)]]))
    print_res(zip(range(3, 4*n, 4), sums))

    print "Getting sums of diagonal links"
    sums = []
    for force in range(5, 4*n - 3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:,f(force)]]))
    print_res(zip(range(5, 4*n, 4), sums))

    for row in xrange(m.shape[0]):
        for column in xrange(m.shape[1]):
            if m[row, column] != 0.0:
                print ("Vertex " + str(row / 2) + " " + ("horizontal" if row % 2 == 0 else "vertical") + ": " +
                    force_position(column))


def parse_input():
    import sys
    filename = sys.argv[1]
    data = open(filename, 'r')
#show_matrix(m[:, 0:np.shape(m)[1] - 1])

#show_matrix(m)
#dim_n = np.shape(m)[0]
#square = m[:, 0:dim_n]  # Take the last row out
#from ipdb import set_trace; set_trace()
#show_matrix(square)
#B = m[:, n]  # Last row
#check_matrix(square)
# fact, perm = lu(square)
# solution = solve_lu(fact, perm, B)
# show_matrix(solution)
#resolved = gauss(square, B)
format_result(solve_problem(m))


# <codecell>


# <codecell>
