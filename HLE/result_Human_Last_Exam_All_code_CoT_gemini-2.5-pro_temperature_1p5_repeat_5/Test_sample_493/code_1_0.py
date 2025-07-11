import math

# The theoretical formula for the average number of stars per constellation (C) is:
# C = 8/3 + sqrt(3)/pi

# Calculate the first term of the equation
term1_numerator = 8
term1_denominator = 3
term1 = term1_numerator / term1_denominator

# Calculate the second term of the equation
# Numerator is the square root of 3
term2_numerator = math.sqrt(3)
# Denominator is pi
term2_denominator = math.pi
term2 = term2_numerator / term2_denominator

# The final result is the sum of the two terms
average_constellation_size = term1 + term2

# Output the equation and the calculated values
print("The average number of stars per constellation (C) is derived from the formula:")
print(f"C = {term1_numerator}/{term1_denominator} + sqrt(3)/pi")
print("\nCalculating the numerical value for each term:")
print(f"Term 1 ({term1_numerator}/{term1_denominator}) = {term1}")
print(f"Term 2 (sqrt(3)/pi) = {term2}")
print("\nSumming the terms gives the final average size:")
print(f"C = {term1} + {term2}")
print(f"\nAverage number of stars per constellation = {average_constellation_size}")
