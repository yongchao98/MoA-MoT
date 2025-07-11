import sympy

# Define 'a' as a symbolic variable representing the positive real constant.
a = sympy.Symbol('a', real=True, positive=True)

# Based on the derivation, the total arc length L for the curve defined by
# the condition 0 <= x <= a is given by the expression L = 3 * a^(2/3).
# L is the sum of the lengths of the symmetric branches in the first and fourth quadrants.
L_expression = 3 * a**(sympy.S(2)/3)

# The problem states that the arc length is 3/2.
given_L = sympy.S(3)/2

# Set up the equation to find the value of 'a'.
equation = sympy.Eq(L_expression, given_L)

# Solve the equation for 'a'.
solution = sympy.solve(equation, a)

# Since 'a' is a positive real number, we take the positive real solution.
# The result of the solve is a list, so we extract the first element.
final_a_value = solution[0]

# Print the final equation with its numerical components.
# The equation is 3 * a^(2/3) = 3/2.
coeff = 3
exp_num = 2
exp_den = 3
rhs_num = 3
rhs_den = 2

print("The final equation relating the constant 'a' to the given arc length is:")
print(f"{coeff} * a^({exp_num}/{exp_den}) = {rhs_num}/{rhs_den}")

# Print the step-by-step solution.
print("\nSolving this equation for 'a':")
print(f"a^({exp_num}/{exp_den}) = ({rhs_num}/{rhs_den}) / {coeff}")
print(f"a^({exp_num}/{exp_den}) = 1/2")
print(f"a = (1/2)^({exp_den}/{exp_num})")

# Print the final exact and numerical values for 'a'.
print(f"\nThe possible value for the constant 'a' is:")
print(f"Exact value: {final_a_value}")
print(f"Numerical value: {final_a_value.evalf()}")
