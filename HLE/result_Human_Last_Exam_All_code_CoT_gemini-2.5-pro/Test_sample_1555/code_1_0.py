import sympy

# The quantum number for the lightest solitonic excitation.
k_lightest = 1

# The quantum number for the subsequent higher excitation.
k_next = 2

# Define N as a symbolic variable for the large-N limit.
N = sympy.Symbol('N')

# In the CP(N-1) model, the mass of the k-th kink-antikink bound state is
# given by M_k proportional to sin(pi*k/N). We form the ratio for the
# excitations corresponding to k_next and k_lightest.
mass_ratio_expression = sympy.sin(k_next * sympy.pi / N) / sympy.sin(k_lightest * sympy.pi / N)

# Compute the limit of the mass ratio as N approaches infinity.
final_ratio = sympy.limit(mass_ratio_expression, N, sympy.oo)

# The final equation for the asymptotic mass ratio is M_2 / M_1 = 2.
# As requested, here are the numbers involved in determining this result.
print(f"Quantum number for the lightest excitation (k_1): {k_lightest}")
print(f"Quantum number for the subsequent excitation (k_2): {k_next}")
print(f"The final asymptotic mass ratio (M_2 / M_1) is: {final_ratio}")