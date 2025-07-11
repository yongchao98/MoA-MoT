import sympy

# The problem asks for the maximum number of vertices labelled 'r' within the interval ]0, 1[.
# A vertex 'v' is part of the dessin d'enfant Gamma if its image phi(v) is in the interval [0, 1].
# A vertex labelled 'r' is defined as a preimage of infinity, so phi(v_r) = infinity.
# It is impossible for a vertex to satisfy both conditions simultaneously: phi(v) cannot be both in [0, 1] and be infinity.
# Therefore, no vertex labelled 'r' can be part of the dessin Gamma.
# The maximum number of such vertices within ]0, 1[ is 0.

max_number_of_r_vertices = 0

# The equation is trivial in this case as it's a direct result of definitions.
# We will print the number '0' as part of a symbolic equation to adhere to the format.
# Let N_r be the number of r-vertices in ]0, 1[.
N_r = sympy.Symbol('N_r')
equation = sympy.Eq(N_r, max_number_of_r_vertices)

# We are asked to output each number in the final equation.
# The equation is N_r = 0.
print(f"{equation.lhs} = {equation.rhs}")
