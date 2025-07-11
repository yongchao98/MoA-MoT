import sympy

# Let m_squared be a symbolic variable representing m^2
m_squared = sympy.Symbol('m_squared')

# The problem states that the 5 degrees of freedom from the tensor part have a squared mass of m^2.
# Let's call this M2_tensor.
M2_tensor = m_squared

# From our analysis of the decomposed Lagrangian, we found the theoretical expression for the squared mass of the tensor modes:
# M2_tensor_theory = -m_squared / C
# where C is an unknown normalization constant of the kinetic Lagrangian.
C = sympy.Symbol('C')
M2_tensor_theory = -m_squared / C

# By equating the given mass with the theoretical expression, we can solve for C.
# m_squared = -m_squared / C
equation_for_C = sympy.Eq(M2_tensor, M2_tensor_theory)
# This simplifies to 1 = -1 / C, so C = -1.
C_val = sympy.solve(equation_for_C, C)[0]


# Our analysis also gave the theoretical expression for the squared mass of the scalar mode (the 6th DoF):
# M2_scalar_theory = m_squared / C
M2_scalar_theory = m_squared / C

# Now substitute the value of C we found back into the expression for the scalar's squared mass.
M2_scalar_final = M2_scalar_theory.subs(C, C_val)

# The result is the squared mass of the sixth degree of freedom.
print("The squared mass of the tensor degrees of freedom (given) is m^2.")
print(f"The analysis of the Lagrangian leads to M_tensor^2 = -m^2 / C. Equating these gives C = {C_val}.")
print(f"The squared mass of the scalar degree of freedom is M_scalar^2 = m^2 / C.")
print(f"Substituting C = {C_val}, the squared mass of the sixth degree of freedom is: {M2_scalar_final}")
# For the final equation format requested
print(f"\nFinal Equation: M_6th^2 = ({sympy.pretty(M2_scalar_final/m_squared)}) * m^2")
