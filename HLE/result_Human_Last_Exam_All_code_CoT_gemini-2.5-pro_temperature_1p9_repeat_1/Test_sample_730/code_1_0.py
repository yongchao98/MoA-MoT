from fractions import Fraction

# The problem is to estimate the mass of a rock with a given radius and density.
# The calculation must use integers no larger than 10.

# Given values as fractions:
# Radius (r) = 0.5 cm = 1/2 cm
# Density (ρ) = 0.9 kg/cm^3 = 9/10 kg/cm^3
r_num, r_den = 1, 2
d_num, d_den = 9, 10

# Volume formula for a sphere contains the fraction 4/3.
v_const_num, v_const_den = 4, 3

# We need a simple fractional approximation for pi that keeps the calculation error under 10%.
# This requires our pi approximation to be between roughly 2.83 and 3.46.
# The simplest integer choice is 3.
pi_approx = 3

# Now, we construct the full mass calculation equation:
# mass = density * volume = (9/10) * (4/3) * pi_approx * (1/2)^3
# The integers involved are {9, 10, 4, 3, 3, 1, 2}. The largest is 10, which is allowed.

# Perform the calculation using the fractions module.
density = Fraction(d_num, d_den)
volume_constant = Fraction(v_const_num, v_const_den)
radius = Fraction(r_num, r_den)
pi_approximation = Fraction(pi_approx)

# The calculation:
mass = density * volume_constant * pi_approximation * (radius**3)

# Print the final equation with each number.
print(f"The parrot's calculation is:")
print(f"mass = ({d_num}/{d_den}) * ({v_const_num}/{v_const_den}) * {pi_approx} * ({r_num}/{r_den})^3")
print(f"The estimated mass is {mass.numerator}/{mass.denominator} kg.")
