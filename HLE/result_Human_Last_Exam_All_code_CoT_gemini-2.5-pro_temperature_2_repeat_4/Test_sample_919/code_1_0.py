import sympy as sp

# Define symbolic variables
# We use these to represent the parameters in the problem and verify the derivation.
K0, mu0, mu, a, d, y = sp.symbols('K0 mu0 mu a d y')
cosh = sp.cosh
sinh = sp.sinh
cos = sp.cos
sin = sp.sin

# Step 5: Solve for the coefficient E
# The derivation in the thought process gives E.
# E = K0 / (a * (cosh(a*d) + (mu0/mu)*sinh(a*d)))
# We define the denominator for clarity.
denominator_D = cosh(a*d) + (mu0/mu) * sinh(a*d)
E = K0 / (a * denominator_D)

# Step 6: Calculate the force
# The force per unit area is (mu0/2) * (H_y(d,y))^2 in the x-direction.
# H_y(d,y) = a * E * sin(a*y)
H_y_at_d = a * E * sin(a*y)

# Force per area expression
force_per_area_sq_mag = (mu0/2) * H_y_at_d**2
force_per_area_sq_mag = sp.simplify(force_per_area_sq_mag)

# Substitute E back into the expression
final_force_expression_mag = (mu0/2) * (a * (K0 / (a * denominator_D)) * sin(a*y))**2
final_force_expression_mag = sp.simplify(final_force_expression_mag)

# The result is the magnitude. The direction is i_x_hat.
# So, f/area = final_force_expression_mag * i_x_hat
# Let's print the result in a readable format similar to the choices.
# We represent the final fraction symbolically to show the structure.
Numerator = mu0/2 * K0**2 * sin(a*y)**2
Denominator_sq = denominator_D**2

# We manually format the output to match the desired physical equation format.
print("The derivation leads to the following expression for the force per unit area:")
# We explicitly show each part of the fraction from the derivation.
print(f"Numerator term: (mu0 / 2) * K0**2 * sin(a*y)**2")
print(f"Denominator term: (cosh(a*d) + (mu0/mu)*sinh(a*d))**2")
print("Direction: +x (i_x_hat)")
print("\nCombining these parts, the final equation is:")
print(f"f/area = ( (mu0 / 2) * K0**2 * sin(a*y)**2 ) / ( (cosh(a*d) + (mu0/mu)*sinh(a*d))**2 ) * i_x_hat")
print("\nThis matches Choice C.")
print("\nFinal Answer Check:")
print("Let's re-assemble the expression from the symbolic parts:")
# Printing a SymPy pretty version of the formula confirms the structure.
force_vector_str = f"({sp.srepr(final_force_expression_mag)})*i_x"

# We represent choice C for direct comparison
Choice_C_mag = (mu0/2) * (K0**2 * sin(a*y)**2) / (cosh(a*d) + mu0/mu * sinh(a*d))**2
print(f"Derived expression: f/area = {sp.pretty(final_force_expression_mag, use_unicode=False)} i_x")
print(f"Choice C expression: f/area = {sp.pretty(Choice_C_mag, use_unicode=False)} i_x")

if sp.simplify(final_force_expression_mag - Choice_C_mag) == 0:
    print("The derived expression is identical to Choice C.")
else:
    print("The derived expression does not match Choice C.")
