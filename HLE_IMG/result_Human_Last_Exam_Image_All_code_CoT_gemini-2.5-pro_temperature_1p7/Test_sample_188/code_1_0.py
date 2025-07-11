import sympy

# Define the symbolic variable 'k' for the number of colors.
k = sympy.Symbol('k')

# Polynomial for a 4-cycle (C4), which is the graph G/(3,4).
P_C4 = (k - 1)**4 + (k - 1)
P_C4_expanded = sympy.expand(P_C4)

# Polynomial for C4 with a pendant vertex, which is the graph G-(3,4)-(2,3).
P_G2 = P_C4 * (k - 1)
P_G2_expanded = sympy.expand(P_G2)

# Polynomial for K4 with one edge removed, which is the graph G-(3,4)/(2,3).
P_H = k * (k - 1) * (k - 2)**2
P_H_expanded = sympy.expand(P_H)

# Polynomial for G-(3,4) by subtracting the two polynomials above.
P_G1 = P_G2_expanded - P_H_expanded
P_G1_expanded = sympy.expand(P_G1)

# Final polynomial for the original graph G.
P_G = P_G1_expanded - P_C4_expanded
P_G_expanded = sympy.expand(P_G)

print("The Chromatic Polynomial of the given graph is P(G, k).")
print("Using the deletion-contraction method, we arrive at the final polynomial.")
print("\nThe final Chromatic Polynomial is:")

# We format the output string to explicitly show each coefficient as requested.
print(f"P(G, k) = {1}*k^5 - {7}*k^4 + {19}*k^3 - {23}*k^2 + {10}*k")

# Also show the factored form for insight.
P_G_factored = sympy.factor(P_G_expanded)
print("\nIn factored form, the polynomial is:")
print(f"P(G, k) = {str(P_G_factored).replace('**', '^')}")
