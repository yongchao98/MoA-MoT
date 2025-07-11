from fractions import Fraction

# The problem requires calculating the flux through two opposite sides of a pyramid.
# We choose the left and right sides. The calculation for one side (the right side)
# leads to evaluating a definite integral which splits into two parts.
# Flux_one_side = Integral_1 + Integral_2
#
# Part 1 of the integral, from the term 2*z*(1-z/4)^6, evaluates to 4/7.
integral_part_1 = Fraction(4, 7)
#
# Part 2 of the integral, from the term z/2 - z^2/8, evaluates to 4/3.
integral_part_2 = Fraction(4, 3)

# The flux through one side is the sum of these two parts.
flux_one_side = integral_part_1 + integral_part_2

# The flux through the other chosen yellow side is identical.
# So, the total flux is twice the flux through one side.
total_flux = 2 * flux_one_side

print("The calculation for the total energy flow is based on the sum of definite integrals.")
print(f"Flux for one yellow side = {integral_part_1.numerator}/{integral_part_1.denominator} + {integral_part_2.numerator}/{integral_part_2.denominator} = {flux_one_side.numerator}/{flux_one_side.denominator}")
print(f"Total flux through two yellow sides = 2 * ({flux_one_side.numerator}/{flux_one_side.denominator}) = {total_flux.numerator}/{total_flux.denominator}")
print(f"The numerical value is approximately {float(total_flux):.4f}")