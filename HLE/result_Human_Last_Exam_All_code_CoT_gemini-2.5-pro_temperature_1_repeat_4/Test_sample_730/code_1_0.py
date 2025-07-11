from fractions import Fraction

# Define the components of the calculation based on the problem and the parrot's constraints.
# Volume formula constant
volume_constant = Fraction(4, 3)
# Radius r = 0.5 cm is represented as a fraction
radius = Fraction(1, 2)
# Density rho = 0.9 kg/cm3 is represented as a fraction
density = Fraction(9, 10)
# Pi is approximated as 3/1 to keep the error below 10% and use integers <= 10.
pi_approximation = Fraction(3, 1)

print("Yes, the parrot can estimate the mass with the given constraints.")
print("Here is the calculation using integers no larger than 10:\n")

# Print the formula with the values substituted
print(f"Mass = (Volume Constant) * (pi_approximation) * (radius)^3 * (density)")
print(f"Mass = ({volume_constant.numerator}/{volume_constant.denominator}) * ({pi_approximation.numerator}/{pi_approximation.denominator}) * ({radius.numerator}/{radius.denominator})^3 * ({density.numerator}/{density.denominator})")

# Expand the cube of the radius
radius_cubed = radius**3
print("\nAfter calculating the cube of the radius:")
print(f"Mass = ({volume_constant.numerator}/{volume_constant.denominator}) * ({pi_approximation.numerator}/{pi_approximation.denominator}) * ({radius_cubed.numerator}/{radius_cubed.denominator}) * ({density.numerator}/{density.denominator})")

# Multiply the numerators and denominators to show the unsimplified fraction
unsimplified_numerator = volume_constant.numerator * pi_approximation.numerator * radius_cubed.numerator * density.numerator
unsimplified_denominator = volume_constant.denominator * pi_approximation.denominator * radius_cubed.denominator * density.denominator
print("\nMultiplying the fractions together:")
print(f"Mass = {unsimplified_numerator} / {unsimplified_denominator}")

# Calculate the final simplified mass
final_mass = volume_constant * pi_approximation * radius_cubed * density
print("\nSimplifying the fraction gives the final estimated mass:")
print(f"Mass = {final_mass.numerator}/{final_mass.denominator} kg")