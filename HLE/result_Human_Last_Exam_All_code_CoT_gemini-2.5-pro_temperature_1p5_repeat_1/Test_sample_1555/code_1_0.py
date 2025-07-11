import sympy

# In the large-N limit of the CP(N-1) model, the mass spectrum of stable
# bound states is given by the well-known formula derived by Witten and others:
# m_k = m_1 * sin(k*pi/N) / sin(pi/N)
# where k is an integer (1, 2, ...), m_1 is the mass of the lightest particle,
# and N is the parameter from the CP(N-1) group.

# The problem asks for the mass ratio between the lightest excitation (k=1)
# and the subsequent higher excitation (k=2) in the limit N -> infinity.

# For k=1, the mass is m_1.
# For k=2, the mass is m_2 = m_1 * sin(2*pi/N) / sin(pi/N).

# The required ratio is m_2 / m_1, which is:
# ratio = sin(2*pi/N) / sin(pi/N)

# We use symbolic math to compute the limit of this ratio as N approaches infinity.

# Define N as a symbol for the limit calculation
N = sympy.Symbol('N')

# Define the expression for the mass ratio
mass_ratio_expression = sympy.sin(2 * sympy.pi / N) / sympy.sin(sympy.pi / N)

# To find the asymptotic value, we calculate the limit of the expression as N -> infinity.
# Using the trigonometric identity sin(2x) = 2*sin(x)*cos(x), the expression simplifies to:
# ratio = 2 * cos(pi/N)
# The limit is lim_{N->inf} 2 * cos(pi/N) = 2 * cos(0) = 2 * 1 = 2.

# We will now use sympy to compute this and format the output as requested.

# Calculate the final asymptotic ratio
asymptotic_ratio = sympy.limit(mass_ratio_expression, N, sympy.oo)

# The limit calculation can be broken down into factors.
# The factor '2' comes from the double-angle formula for sine.
# The other factor is the limit of cos(pi/N) as N -> infinity.
factor1 = 2
factor2 = sympy.limit(sympy.cos(sympy.pi / N), N, sympy.oo)

print("The asymptotic mass ratio m_2/m_1 is found by taking the limit of the expression as N approaches infinity.")
print("The final calculation step is based on the simplified form 2 * cos(\u03C0/N):")
print(f"{factor1} * {int(factor2)} = {int(asymptotic_ratio)}")