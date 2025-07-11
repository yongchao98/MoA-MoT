import sympy

# Define the variable k for our polynomial
k = sympy.Symbol('k')

# Method 1: Deletion-Contraction using e = (1,5)
# Chromatic polynomial of the graph with edge (1,5) deleted
P_G_minus_e = k * (k - 1)**2 * (k - 2)**2

# Chromatic polynomial of the graph with edge (1,5) contracted (which is K4)
P_G_contract_e = k * (k - 1) * (k - 2) * (k - 3)

# The chromatic polynomial of the original graph G
P_G = P_G_minus_e - P_G_contract_e

# Factor the resulting polynomial to get the simplified form
factored_poly = sympy.factor(P_G)

# We can also verify this using a direct counting argument which gives:
# P_G_direct = k*(k-1)*(k-2) * ( (k-1) + (k-2)*(k-3) )
# factored_poly_direct = sympy.factor(P_G_direct)
# assert factored_poly == factored_poly_direct

# The final polynomial in factored form
# k * (k-1) * (k-2) * (k**2 - 4*k + 5)
# Let's format this for a clear output
final_poly_str = "k * (k - 1) * (k - 2) * (k**2 - 4*k + 5)"

print("The chromatic polynomial for the given graph is:")
print(f"P(k) = {final_poly_str}")

# The polynomial can also be expanded:
# expanded_poly = sympy.expand(factored_poly)
# print("\nExpanded form:")
# print(f"P(k) = {expanded_poly}")
# Output: P(k) = k**5 - 7*k**4 + 19*k**3 - 23*k**2 + 10*k
