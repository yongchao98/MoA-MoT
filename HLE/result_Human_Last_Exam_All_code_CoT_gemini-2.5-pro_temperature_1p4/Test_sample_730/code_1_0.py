from fractions import Fraction

# Given values and constants represented as integers
density_num = 9
density_den = 10
radius_num = 1
radius_den = 2
vol_const_num = 4
vol_const_den = 3
exponent = 3

# The parrot's approximation for pi (A/B, where A, B <= 10)
# We use 3/1 as it satisfies the error constraint (<10%)
pi_approx_num = 3
pi_approx_den = 1

# Calculate the final mass as a fraction
mass_approx = (Fraction(density_num, density_den) *
               Fraction(vol_const_num, vol_const_den) *
               Fraction(pi_approx_num, pi_approx_den) *
               (Fraction(radius_num, radius_den) ** exponent))

# Print the final equation with all numbers involved
radius_cubed_num = radius_num ** exponent
radius_cubed_den = radius_den ** exponent
print("The parrot's calculation is:")
print(f"Mass = ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({radius_cubed_num}/{radius_cubed_den})")
print(f"Mass = {mass_approx.numerator}/{mass_approx.denominator} kg")