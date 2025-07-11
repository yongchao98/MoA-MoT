import sympy

# The purpose of this script is to present the final solution to the differential equation.
# The solving process involves advanced techniques (transformation to a Lagrange equation)
# which are explained in the text description. This code simply formats and prints the result.

# The general solution to the differential equation is y^2 = C*x^2 + C^2 + C + 9,
# where C is an arbitrary constant.

# We define the numbers present in the final equation as variables to explicitly show them.
power_of_y = 2
power_of_x = 2
power_of_C = 2
coefficient_of_C = 1
constant_term = 9

# We use sympy to pretty print the equation, but it's not necessary for the result itself.
# Let's build the equation string manually to match the required output format.
# "y**2 = C * x**2 + C**2 + C + 9"

print("The general solution of the equation is:")
print(f"y**{power_of_y} = C * x**{power_of_x} + C**{power_of_C} + {coefficient_of_C}*C + {constant_term}")
print("\nWhere C is an arbitrary constant.")
