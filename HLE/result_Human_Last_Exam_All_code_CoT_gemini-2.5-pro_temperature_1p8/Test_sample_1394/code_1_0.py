import sympy

# Define the variables
y, x = sympy.symbols('y x')
C = sympy.Symbol('C')  # Arbitrary constant

# The general solution equation
# y**2 = C*x**2 + C**2 + C + 9
# Rearrange to y**2 - C*x**2 - C**2 - C - 9 = 0
equation = y**2 - C*x**2 - C**2 - C - 9

# Print the solution in a human-readable format
print("The general solution to the differential equation is:")
solution_str = f"y^2 = C*x^2 + C^2 + C + 9"
print(solution_str)
print("where C is an arbitrary constant.")
print("\nHere are the numbers present in the final equation (y^2 - C*x^2 - C^2 - C - 9 = 0):")

# Extract and print the numbers (powers and constants)
# For y^2
power_y = 2
print(f"The power of y is: {power_y}")

# For C*x^2
power_x = 2
print(f"The power of x is: {power_x}")

# For C^2
power_C_1 = 2
print(f"The power of C in the 'C^2' term is: {power_C_1}")

# For C
power_C_2 = 1
print(f"The power of C in the 'C' term is: {power_C_2}")

# For the constant term 9
constant_term = 9
print(f"The constant term is: {constant_term}")
