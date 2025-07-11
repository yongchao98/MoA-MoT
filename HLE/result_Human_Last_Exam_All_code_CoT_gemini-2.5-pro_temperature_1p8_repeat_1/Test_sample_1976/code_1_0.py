import sympy

# Use sympy for symbolic manipulation and better printing
n = sympy.Symbol('n', integer=True, odd=True)

# Step 1: Normalization
# Tr(phi_2) = 1
# Tr(phi_2^perp) = Tr(I_4 - phi_2) = 4 - 1 = 3
# Z_n = Tr(bigotimes(phi_2)) + (1/3) * Tr(bigotimes(phi_2^perp))
#     = Tr(phi_2)^(n+1) + (1/3) * Tr(phi_2^perp)^(n+1)
#     = 1^(n+1) + (1/3) * 3^(n+1) = 1 + 3^n
Z_n = 1 + 3**n

print("Step 1: The normalization constant Z_n for J_n is Tr(J_n_unnormalized).")
print(f"Z_n = 1 + 3^n = {sympy.pretty(Z_n)}\n")

# Step 2 & 3: Correlation matrix entries
# Let's use the Pauli basis sigma_a = sigma_{a_1} x ... x sigma_{a_{n+1}}
# t_ab = Tr(J_n * (sigma_a tensor sigma_b))
# The state and operators factorize over the n+1 pairs.
# t_ab = (1/Z_n) * [ Prod_k Tr(phi_2 * (s_{a_k} x s_{b_k})) + (1/3) * Prod_k Tr(phi_2^perp * (s_{a_k} x s_{b_k})) ]
# Let C_{ij} = Tr(phi_2 * (s_i x s_j)) and D_{ij} = Tr(phi_2^perp * (s_i x s_j)).
# C is diagonal: C = diag(1, 1, -1, 1) for i,j in {0,1,2,3}
# D is diagonal: D = diag(3, -1, 1, -1)
# Therefore, t_ab is non-zero only if a_k = b_k for all k, so T is diagonal.
# We need to compute t_aa for a != 0 vector.
print("Step 2 & 3: The correlation matrix T is diagonal in the Pauli basis.")
print("The diagonal entries t_aa for a = (a_1, ..., a_{n+1}) are calculated.\n")

# Step 4: Sum the absolute values
# We sum over all a != 0. Let n_i be the number of times i appears in a.
# t_aa = (1/Z_n) * [ (-1)^n_2 + (1/3) * 3^n_0 * (-1)^(n_1+n_3) ]
# Since n is odd, we have two main cases for the term inside |...|:
# 1. n_0 is even: n_1+n_2+n_3 = n+1-n_0 is even. Then n_2 and n_1+n_3 have the same parity.
#    |s + 3^(n_0-1)s| = 1 + 3^(n_0-1) for n_0>=2 and 4/3 for n_0=0.
# 2. n_0 is odd: n_1+n_2+n_3 = n+1-n_0 is odd. n_2 and n_1+n_3 have opposite parity.
#    |s - 3^(n_0-1)s| = 3^(n_0-1) - 1 for n_0>=3 and 0 for n_0=1.

# Let's compute the total sum of absolute values, denoted Sigma.
# We sum S_n0 over n_0 from 0 to n.
# S_n0 = C(n+1, n_0) * 3^(n+1-n_0) * V(n_0), where V is the value above.
print("Step 4: Summing the absolute values of the diagonal entries.")

# For a specific odd n to demonstrate
n_val = 3
print(f"Let's demonstrate for n = {n_val}:\n")

Z_val = 1 + 3**n_val
total_sum = 0
for n0 in range(n_val + 1):
    # Number of ways to choose positions for n_0 identities
    num_configs = sympy.binomial(n_val + 1, n0)
    # Number of ways to fill the rest with sigma_x,y,z
    num_pauli_strings = 3**(n_val + 1 - n0)
    num_vectors = num_configs * num_pauli_strings

    value = 0
    if n0 % 2 == 1: # n0 is odd
        if n0 == 1:
            value = 0
        else: # n0 >= 3
            value = 3**(n0 - 1) - 1
    else: # n0 is even
        if n0 == 0:
            value = sympy.Rational(4, 3)
        else: # n0 >= 2
            value = 1 + 3**(n0 - 1)

    term_sum = num_vectors * value
    total_sum += term_sum
    print(f"For n0 = {n0}: sum term = C({n_val+1},{n0}) * 3^({n_val+1-n0}) * {value} = {term_sum}")

print(f"\nTotal sum Sigma = {total_sum}")
norm = total_sum / Z_val
print(f"Normalization Z_n = {Z_val}")
print(f"The 1-norm ||T||_1 = Sigma / Z_n = {total_sum} / {Z_val} = {norm}\n")


# Step 5: Final Result
# The sum can be simplified for a general odd n to:
# Sigma = (2**(n+1) - 1) * (3**n + 1)
# So, ||T||_1 = Sigma / Z_n = ( (2**(n+1)-1)*(3**n+1) ) / (3**n+1) = 2**(n+1)-1
n = sympy.Symbol('n')
final_norm = 2**(n+1) - 1

print("Step 5: For a general odd n, the sum simplifies beautifully.")
print(f"The final expression for the 1-norm is: {sympy.pretty(final_norm)}")
final_eq_lhs = sympy.Symbol("||T||_1")
final_eq = sympy.Eq(final_eq_lhs, final_norm)
n_val = 3
print(f"For n = {n_val}: ||T||_1 = 2^({n_val}+1) - 1 = {2**(n_val+1)-1}")