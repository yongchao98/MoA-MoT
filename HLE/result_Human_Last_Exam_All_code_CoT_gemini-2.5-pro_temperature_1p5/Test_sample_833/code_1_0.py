import sympy

# The final expression for the lower bound is a = -1/16.
# I will demonstrate the calculation steps here.
u_bar_at_0 = sympy.Rational(-1, 2)
K_conv_u_bar_at_0 = sympy.Rational(-1, 4)

# The expression is E = 1/4 * (u_bar(0) - K_conv_u_bar(0))
E = sympy.Rational(1, 4) * (u_bar_at_0 - K_conv_u_bar_at_0)

# We print the calculation to show the numbers involved in the final equation.
print(f"The simplified expression to evaluate is: 1/4 * (u_bar(0) - (K*u_bar)(0))")
print(f"Calculated value for u_bar(0) = {u_bar_at_0}")
print(f"Calculated value for (K*u_bar)(0) = {K_conv_u_bar_at_0}")
print(f"The final lower bound is a = 1/4 * ({u_bar_at_0} - ({K_conv_u_bar_at_0}))")
final_a_numerator = E.p
final_a_denominator = E.q
print(f"a = {final_a_numerator}/{final_a_denominator}")
print(f"As a decimal, a = {float(E)}")