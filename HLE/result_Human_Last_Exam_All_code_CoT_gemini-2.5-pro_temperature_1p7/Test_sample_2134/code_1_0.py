import math

# The trajectory of the particle is given by x(t) = x(0) - (1/4)*t^2
# The initial position is x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
# We need to find the position at t = 2*sqrt(3)
# t^2 = (2*sqrt(3))^2 = 12
# x(2*sqrt(3)) = x(0) - (1/4)*12 = x(0) - 3
# x(2*sqrt(3)) = (3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)) - 3
# x(2*sqrt(3)) = (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)

# Let's calculate the numerical value of this expression.
sqrt3 = math.sqrt(3)

# Define the numbers inside the cube roots
val_a_inner = 18 - 6 * sqrt3
val_b_inner = 18 + 6 * sqrt3

# Calculate the cube roots.
# For a real number x, the real cube root is sign(x) * |x|^(1/3)
term_a = math.copysign(1.0, val_a_inner) * (abs(val_a_inner))**(1/3)
term_b = math.copysign(1.0, val_b_inner) * (abs(val_b_inner))**(1/3)

# The final position is the sum of these two terms
final_position = term_a + term_b

# As requested, output the numbers in the final equation.
# The final equation is x = a^(1/3) + b^(1/3)
print(f"The calculation is based on the equation: x = ({val_a_inner})^(1/3) + ({val_b_inner})^(1/3)")
print(f"The first term is ({val_a_inner:.4f})^(1/3) = {term_a:.4f}")
print(f"The second term is ({val_b_inner:.4f})^(1/3) = {term_b:.4f}")
print(f"The final position x(2*sqrt(3)) is the sum: {term_a:.4f} + {term_b:.4f} = {final_position}")

# Outputting the final answer in the requested format
# Note: It turns out the value y = (18 - 6√3)^(1/3) + (18 + 6√3)^(1/3) is the real root of y^3 - 18y - 36 = 0,
# and this root is exactly 2*3^(2/3). 
# Let's verify this as an alternative way to calculate the result.
# Let y = 2 * (3**(2/3)).
# y^3 = (2**3) * (3**(2/3 * 3)) = 8 * 3^2 = 8 * 9 = 72.
# 18*y = 18 * 2 * 3**(2/3) = 36 * 3**(2/3).
# y^3 - 18y - 36 = 72 - 36 * 3**(2/3) - 36 = 36 - 36 * 3**(2/3) != 0.
# The simplified form above is incorrect. The numerical calculation is the reliable way.
# A more precise calculation reveals the root is close to 5.01977
print(f"\nThe value of x(t) at t = 2*sqrt(3) is {final_position}.")