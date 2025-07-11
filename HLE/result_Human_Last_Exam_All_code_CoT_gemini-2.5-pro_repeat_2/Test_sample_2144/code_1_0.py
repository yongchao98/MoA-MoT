import math

# The value of the parameter p where y = -3 was found to be 3.
p_star = 3

# The constant of integration was found to be 3.
C = 3

# The equation for the x-coordinate of the trajectory is x(p) = -(3p + C) / sqrt(2p + C).
# We substitute p = p_star to find x0.
numerator = -(3 * p_star + C)
denominator = math.sqrt(2 * p_star + C)
x0 = numerator / denominator

print(f"The parameter p where y = -3 is p_star = {p_star}.")
print(f"The position x0 is calculated using the equation: x0 = -(3 * p_star + C) / sqrt(2 * p_star + C)")
print(f"Substituting the values: x0 = -(3 * {p_star} + {C}) / sqrt(2 * {p_star} + {C})")
print(f"x0 = {numerator} / sqrt({2 * p_star + C})")
print(f"x0 = {numerator} / {denominator}")
print(f"The position x0 is: {x0}")
