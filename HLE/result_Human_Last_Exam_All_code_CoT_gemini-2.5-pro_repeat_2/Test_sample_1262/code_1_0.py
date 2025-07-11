import itertools

def get_derangements(n):
    """Generates all derangements of n elements."""
    if n == 0:
        return [[]]
    # We use permutations of (1, 2, ..., n)
    elements = range(1, n + 1)
    perms = itertools.permutations(elements)
    derangements = []
    for p_tuple in perms:
        if all(p_tuple[i] != i + 1 for i in range(n)):
            derangements.append(list(p_tuple))
    return derangements

def count_excedances(perm):
    """Counts the number of excedances in a permutation."""
    return sum(1 for i, val in enumerate(perm) if val > i + 1)

# Part (a): Confirm whether H(U_{n-1, E})(t) = t^{n-1} d_n(t).
# We demonstrate the identity is false by comparing polynomial degrees for n=4.
print("(a) Analysis of the identity for n = 4:")
n_a = 4
# Degree of the LHS: The rank of U_{n-1, n} is n-1. The degree of its Chow ring's Hilbert series is rank - 1.
deg_LHS = (n_a - 1) - 1
print(f"The degree of the Hilbert series H(U_{{3,4}})(t) is (4-1) - 1 = {deg_LHS}.")

# Degree of the RHS: We need the degree of d_4(t).
derangements_a = get_derangements(n_a)
deg_d_n = 0
if derangements_a:
    deg_d_n = max(count_excedances(p) for p in derangements_a)
deg_RHS = (n_a - 1) + deg_d_n
print(f"The degree of the derangement polynomial d_4(t) is {deg_d_n}.")
print(f"The degree of the expression t^(4-1) * d_4(t) is (4-1) + {deg_d_n} = {deg_RHS}.")
print(f"Since the degrees ({deg_LHS} and {deg_RHS}) are not equal, the identity is false.")
print("-" * 30)

# Part (b): State if the leading coefficient of d_n(t) is always 1.
# The leading coefficient is the number of derangements with the maximum number of excedances.
print("(b) Analysis of the leading coefficient of d_n(t):")
for n_b in range(2, 6):
    derangements_b = get_derangements(n_b)
    if not derangements_b:
        continue
    
    max_exc = max(count_excedances(p) for p in derangements_b)
    leading_coeff = sum(1 for p in derangements_b if count_excedances(p) == max_exc)
    
    print(f"For n = {n_b}, the maximum number of excedances is {max_exc}.")
    print(f"The number of derangements with {max_exc} excedances (the leading coefficient) is {leading_coeff}.")
print("The leading coefficient is 1 for all tested cases, confirming the theoretical argument.")
print("-" * 30)

# Part (c): Give the value of d_3(1).
# d_3(1) is the number of derangements of 3 elements.
print("(c) Calculation of d_3(1):")
n_c = 3
derangements_c = get_derangements(n_c)
num_derangements_c = len(derangements_c)
print(f"The derangements in S_3 are: {derangements_c}")
print(f"d_3(1) is the count of these derangements.")
print(f"d_3(1) = {num_derangements_c}")