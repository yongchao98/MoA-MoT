# Based on the analytical derivation, the value of the function l(b) is a constant.
# l(b) = 101 for all b in (-1, 1).

# Assign the derived values
l_half = 101
l_neg_half = 101
factor = 6

# Calculate the final result
result = factor * (l_half + l_neg_half)

# Print the final equation with each number
print(f"{factor} * ({l_half} + {l_neg_half}) = {result}")

# Final Answer
print(f"<<<{result}>>>")