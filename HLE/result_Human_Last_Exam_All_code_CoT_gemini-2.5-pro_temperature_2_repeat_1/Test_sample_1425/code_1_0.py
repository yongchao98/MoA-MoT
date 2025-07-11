import sympy as sp

# Define the symbols to be used in the equations
beta, mu, g = sp.symbols('β μ g', real=True)
# Assume beta (inverse temperature) and g (degeneracy) are positive
# beta is 1/(k_B T) and must be positive.
# g is a count of states and must be positive.

# Introduction explaining the setup
print("This script calculates the partition function Z for a system with Hamiltonian H = -μN.")
print("The calculation is based on the path integral formalism in the grand canonical ensemble.")
print("This corresponds to a system with zero single-particle energy and chemical potential μ.")
print("The result depends on the statistics of the constituent particles (bosons or fermions).")
print("Let 'g' be the number of degenerate single-particle states.\n")

print("="*60)
# --- Fermionic Case ---
print("Case 1: The system consists of FERMIONS")
print("-" * 60)
# Define the expression for the fermionic partition function
# The formula is Z = (1 + exp(βμ))^g
Z_fermions_expr = (1 + sp.exp(beta * mu))**g
print("The path integral for fermions evaluates to the following partition function Z_fermions:")
# The sp.pretty() function prints the expression in a formatted, readable way
print(sp.pretty(Z_fermions_expr, use_unicode=True))
print("\n")


print("="*60)
# --- Bosonic Case ---
print("Case 2: The system consists of BOSONS")
print("-" * 60)
# Define the expression for the bosonic partition function
# The formula is Z = (1 - exp(βμ))^(-g)
Z_bosons_expr = (1 - sp.exp(beta * mu))**(-g)
print("The path integral for bosons evaluates to the following partition function Z_bosons:")
print(sp.pretty(Z_bosons_expr, use_unicode=True))
print("\nThis result is valid for μ < 0 to ensure convergence.")
print("="*60)
