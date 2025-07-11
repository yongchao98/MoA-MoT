# Plate parameters
chars_per_circle = 3333
chars_per_square = 90

# Based on the analysis, the optimal strategy is to use a dense packing
# for the high-value circular plates. This allows fitting 9 circles.
# This packing is very tight, leaving no recoverable space for 10x10 squares.

# Optimal number of circular plates
N = 9
# Optimal number of squared plates
M = 0

# Calculate the maximal number of characters
K = N * chars_per_circle + M * chars_per_square

# Output the final result in the specified formats
print("Optimal Configuration (N:M:K):")
print(f"{N}:{M}:{K}")
print("\nCalculation for the maximal number of characters:")
print(f"{N} circular plates * {chars_per_circle} chars/plate + {M} squared plates * {chars_per_square} chars/plate = {K} characters")

# The final answer in the requested format N:M:K
print(f"\n<<<9:0:29997>>>")