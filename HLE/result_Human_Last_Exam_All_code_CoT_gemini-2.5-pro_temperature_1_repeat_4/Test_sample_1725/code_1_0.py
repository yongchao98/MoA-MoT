import math

# This script solves for two quantities related to a continuum model of a crystal.
# The reasoning is laid out in the comments.

# Part A: What is the dimension of pi's fibers?

# The field C is a vector bundle connection on TR^3.
# In a 3-dimensional space, a connection is locally specified by its
# coefficients, C^k_{ij}.
# The indices i, j, and k each range from 1 to 3.
# The total number of independent components is the product of the range of each index.
dim_base_space = 3
num_components = dim_base_space * dim_base_space * dim_base_space

# The fiber of the bundle pi is the space of all possible values of C at a single point.
# Therefore, its dimension is the number of components of C.
dimension_of_fiber = num_components

# Part B: How many coefficients specify E?

# E is a quadratic energy functional of the torsion tensor T^k_{ij} = C^k_{ij} - C^k_{ji}.
# The number of coefficients is the number of independent quadratic invariants of T
# under the cubic symmetry group O_h.
# The torsion tensor T has 9 independent components and transforms under a specific
# 9D representation of O_h, which can be identified as T_1u tensor T_1g.
# The number of quadratic invariants can be found using character theory.
# The character of the symmetric square of the representation for T is computed,
# and its projection onto the trivial representation A_1g gives the number of invariants.
#
# Let V be the 9D representation for the torsion tensor T.
# The number of invariants is n = <Sym^2(V), A_1g>.
# The character for V = T_1u tensor T_1g is (for classes E, 8C3, 6C2, 6C4, 3C2', i, 8S6, 6sigma_d, 6S4, 3sigma_h):
# chi_V = (9, 0, 1, 1, 1, -9, 0, -1, -1, -1)
# The character for V(g^2) is:
# chi_V_g2 = (9, 0, 9, 1, 9, 9, 0, 9, 1, 9)
# The character for the symmetric square Sym^2(V) is chi_S2V = 0.5 * (chi_V^2 + chi_V_g2):
# chi_S2V = (45, 0, 5, 1, 5, 45, 0, 5, 1, 5)
# The number of invariants is the dot product of this character with the A_1g character (all 1s),
# divided by the order of the group O_h (which is 48).
#
# n = (1/48) * [1*45 + 8*0 + 6*5 + 6*1 + 3*5 + 1*45 + 8*0 + 6*5 + 6*1 + 3*5]
# n = (1/48) * [45 + 30 + 6 + 15 + 45 + 30 + 6 + 15]
# n = (1/48) * 192 = 4
#
# The result of this established group-theoretical calculation is 4.
number_of_coefficients = 4

# The final equation is not really an equation, but a printout of the two numbers.
# I will print each number as requested.
print(f"{dimension_of_fiber} {number_of_coefficients}")