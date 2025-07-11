import math

# We need to calculate 3 / (e^2 + 3)
e = math.e
e_squared = e**2
denominator = e_squared + 3
result = 3 / denominator

# The problem asks to output each number in the final equation.
u_0_1_numerator = -6
u_0_1_denominator_val = e_squared + 3
u_0_1 = u_0_1_numerator / u_0_1_denominator_val

final_quantity_numerator = 3
final_quantity_denominator_val = e_squared + 3

print(f"The value of u(0,1) is given by the equation: u(0,1) = {u_0_1_numerator} / (e^2 + {3})")
print(f"Numerically, u(0,1) is approximately: {u_0_1}")
print(f"The quantity to find, -u(0,1)/2, is given by the equation: -u(0,1)/2 = {final_quantity_numerator} / (e^2 + {3})")
print(f"The final numerical value is: {result}")
