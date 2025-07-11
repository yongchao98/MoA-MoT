import sympy

# Plan:
# 1. Define symbolic variables for the coupling (u), the expansion parameter (epsilon),
#    and the series coefficients (B, C).
# 2. Define the one-loop beta function and the one-loop anomalous dimension gamma_m.
# 3. Solve for the non-trivial fixed point u_star.
# 4. Substitute u_star into the expression for nu_inverse = 2 - gamma_m.
# 5. Calculate nu and expand it in epsilon to show the correction term.
# 6. Conclude by stating the order in 'u' from which this correction originates.

# 1. Define symbolic variables
u, epsilon, B, C = sympy.symbols('u epsilon B C', real=True, positive=True)

# 2. Define the one-loop RG functions
# Beta function: beta(u) = -epsilon*u + B*u^2
beta_u = -epsilon * u + B * u**2

# Anomalous dimension for the mass term: gamma_m(u) = C*u
gamma_m_u = C * u

print("Step 1: Define the RG functions (to one loop).")
print(f"Beta function, beta(u) = {beta_u}")
print(f"Anomalous dimension, gamma_m(u) = {gamma_m_u}")
print("-" * 40)

# 3. Solve for the non-trivial fixed point u_star
# We solve beta(u) = 0 for u != 0
u_star_solutions = sympy.solve(beta_u, u)
# The solutions are [0, epsilon/B]. We need the non-trivial one.
u_star = u_star_solutions[1]

print("Step 2: Find the non-trivial fixed point u_star by solving beta(u) = 0.")
print(f"The fixed point is u_star = {u_star}")
print("-" * 40)

# 4. Define the expression for nu_inverse and substitute u_star
# The relation is nu_inverse = 2 - gamma_m(u)
nu_inverse_expr = 2 - gamma_m_u

# Evaluate at the fixed point
nu_inverse_at_u_star = nu_inverse_expr.subs(u, u_star)

print("Step 3: Evaluate nu^-1 = 2 - gamma_m(u) at the fixed point u_star.")
print(f"nu^-1 = {nu_inverse_at_u_star}")
print("-" * 40)

# 5. Calculate nu and its series expansion in epsilon
nu_at_u_star = 1 / nu_inverse_at_u_star

# Expand nu in a series around epsilon = 0 to see the structure
nu_series = nu_at_u_star.series(epsilon, 0, 2)

print("Step 4: Calculate the exponent nu and its series expansion in epsilon.")
print(f"The full expression for nu is: nu = {nu_at_u_star}")
print(f"The series expansion for nu is: nu = {nu_series}")
print("\nThis expansion has the form: 1/2 + (C/(4*B))*epsilon + O(epsilon^2)")
print("The numbers appearing in the final equation are 1, 2, and 4.")
print("-" * 40)

# 6. Explain the result and provide the final answer
print("Conclusion:")
print("The mean-field value for nu (the zeroth-order term) is 1/2.")
print("The first correction to this value comes from the term 'gamma_m(u)' in the expression for nu^-1.")
print(f"The lowest order term in the expansion of gamma_m(u) is '{gamma_m_u}', which is of the first order in the coupling constant 'u'.")

final_order = 1
print(f"\nTherefore, the initial non-vanishing contribution to nu is acquired at order {final_order} in u.")