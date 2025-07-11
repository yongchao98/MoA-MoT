import numpy as np

# The abelianization of the group G, denoted Ab(G), is determined by the
# images of its generators under a homomorphism into Z^2. These images,
# corresponding to the standard generators x_0 and x_1, are the vectors v1 and v2.
v1 = np.array([1, -1])
v2 = np.array([0, -1])

# To find the rank of Ab(G), we can find the rank of the subgroup generated
# by these vectors. We form a matrix M with these vectors as columns.
M = np.array([v1, v2]).T

# The elements of matrix M are:
# M = [[a, b],
#      [c, d]]
a, b = M[0, 0], M[0, 1]
c, d = M[1, 0], M[1, 1]

# The determinant of M is given by the equation: det = ad - bc.
determinant = a * d - b * c

# If the determinant is +/- 1, the vectors generate Z^2, which has rank 2.
# In general, the rank r is the rank of the matrix M.
r = np.linalg.matrix_rank(M)

# The group Ab(G) is isomorphic to Z^r. A free abelian group like Z^r
# has a trivial torsion subgroup {0}. The order of this subgroup, t, is 1.
t = 1

# The problem asks for the pair (r, t). The following prints the values
# and shows the determinant calculation as requested.
print("The determinant of the matrix of generator images is calculated.")
# The following line outputs each number in the final equation for the determinant.
print(f"det = ({a})*({d}) - ({b})*({c}) = {int(determinant)}")
print(f"Since the determinant is non-zero, the vectors are linearly independent.")
print(f"The rank of Ab(G) is r = {int(r)}.")
print(f"The order of the torsion subgroup of Ab(G) is t = {t}.")
