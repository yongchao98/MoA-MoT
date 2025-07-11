import math

# The problem is to find the maximum distance R_max for a volume V=1.
# The relationship derived from maximizing the gravitational field is:
# V = (4 * pi / 15) * R_max^3
# For V = 1, we solve for R_max:
# R_max^3 = 15 / (4 * pi)
# R_max = (15 / (4 * pi))^(1/3)

# Define the numbers in the final equation
numerator = 15
denominator_part1 = 4
denominator_part2 = math.pi
exponent = 1/3

# Print the equation being solved
print("The equation to find the furthest distance (R_max) is:")
print(f"R_max = ({numerator} / ({denominator_part1} * {denominator_part2})) ^ ({exponent})")
print("")

# Calculate the result
r_max_cubed = numerator / (denominator_part1 * denominator_part2)
r_max = r_max_cubed ** exponent

# Print the final answer
print("The furthest point on the surface of the playdough is:")
print(f"{r_max} meters from point A.")