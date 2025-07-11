# Based on the derivation, the value of l(b) is a constant, 101, for any b in (-1, 1).
# We can therefore directly calculate the final expression.

l_val_half = 101
l_val_neg_half = 101
coefficient = 6

# Calculate the final result
result = coefficient * (l_val_half + l_val_neg_half)

# Print the equation with all the numbers
print(f"{coefficient} * ({l_val_half} + {l_val_neg_half}) = {result}")

# The final answer in the requested format
# <<<1212>>>