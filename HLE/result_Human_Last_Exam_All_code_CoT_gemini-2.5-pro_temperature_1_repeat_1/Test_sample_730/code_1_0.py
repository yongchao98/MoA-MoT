from fractions import Fraction

# The problem is to estimate the mass of a sphere with a given radius and density.
# The formula is: Mass = Density * Volume, where Volume = (4/3) * pi * r^3.
# The exact values are: radius r = 0.5 cm = 1/2 cm, density rho = 0.9 kg/cm^3 = 9/10 kg/cm^3.
# The calculation must use integers no larger than 10.

# The true mass is (9/10) * (4/3) * pi * (1/2)^3 = (3 * pi) / 20 kg.
# Using pi ~ 3.14159, the true mass is ~0.4712 kg.
# A 10% error means the estimate must be between ~0.4241 kg and ~0.5183 kg.

# To use the smallest possible integers, we approximate the inputs.
# The parrot can use r = 1/2 and the volume constant 4/3, as these integers (1, 2, 3, 4) are small.
# To minimize the largest integer, we approximate density rho = 0.9 (9/10) with 8/9. This makes the largest integer 9 instead of 10.
# We approximate pi with the whole number 3.

# Let's define these approximations as fractions.
density_approx = Fraction(8, 9)
volume_constant = Fraction(4, 3)
pi_approx = Fraction(3, 1)
radius = Fraction(1, 2)

# The calculation for r^3 is radius * radius * radius
r_cubed_1 = radius
r_cubed_2 = radius
r_cubed_3 = radius

# Now, we calculate the estimated mass using the parrot's friendly numbers.
estimated_mass = density_approx * volume_constant * pi_approx * r_cubed_1 * r_cubed_2 * r_cubed_3

# The estimated mass is 4/9, which is ~0.4444 kg. This is within the 10% error margin of the true mass.
# The integers used are 8, 9, 4, 3, 3, 1, 2. The largest is 9.
# Therefore, the answer is Yes, and the largest integer is 9.

print("Yes, the parrot can perform the estimation. Here is the calculation using integers of 9 or less:")
print(
    f"{density_approx.numerator} / {density_approx.denominator} * "
    f"{volume_constant.numerator} / {volume_constant.denominator} * "
    f"{pi_approx.numerator} * "
    f"{r_cubed_1.numerator} / {r_cubed_1.denominator} * "
    f"{r_cubed_2.numerator} / {r_cubed_2.denominator} * "
    f"{r_cubed_3.numerator} / {r_cubed_3.denominator} = "
    f"{estimated_mass.numerator} / {estimated_mass.denominator}"
)
<<<Y9>>>