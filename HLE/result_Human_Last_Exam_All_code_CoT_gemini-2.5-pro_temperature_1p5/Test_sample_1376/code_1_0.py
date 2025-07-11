# The optimal number of squares (N) and circles (M) has been determined
# through packing analysis.
# N = Number of 10x10cm squares
# M = Number of 20cm radius circles
optimal_N = 0
optimal_M = 9

# Define the number of characters per artifact
chars_per_square = 4
chars_per_circle = 999

# Calculate the total maximum number of characters (K)
# K = (characters from squares) + (characters from circles)
total_chars_K = (chars_per_square * optimal_N) + (chars_per_circle * optimal_M)

# Print the final result and the equation used to calculate it
print("The optimal number of squares (N) is 0.")
print("The optimal number of circles (M) is 9.")
print("\nThe final calculation is:")
# Output each number in the final equation as requested
print(f"{chars_per_square} * {optimal_N} + {chars_per_circle} * {optimal_M} = {total_chars_K}")
print("\nThe result in the format N:M:K is:")
print(f"{optimal_N}:{optimal_M}:{total_chars_K}")