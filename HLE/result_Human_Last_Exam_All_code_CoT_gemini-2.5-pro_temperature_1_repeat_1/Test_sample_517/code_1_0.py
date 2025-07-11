import sympy

# Define symbols
a_01 = sympy.Symbol('a(0,1)')
pi = sympy.pi

# Step 1: The limit p is given by the formula 1 - (2/pi) / a(0,1)
# p = 1 - (2/pi) / a_01
# The value of the potential kernel for a neighbor of the origin is a(0,1) = 4/pi
a_01_val = 4/pi

# Step 2: Substitute the value of a(0,1) into the expression for p
p_val = 1 - (2/pi) / a_01_val

# Step 3: Print the equation and the final result
eq_str = f"p = 1 - (2/\u03C0) / a(0,1) = 1 - (2/\u03C0) / (4/\u03C0) = 1 - 2/4 = {p_val}"
print("The limit is calculated as follows:")
print(eq_str)

final_value = float(p_val)
print(f"\nThe numerical value of the limit is: {final_value}")
