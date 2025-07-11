# Define the given transition rates
l01 = 0.019
l10 = 0.65
l12 = 0.4
l21 = 0.392
l23 = 0.008
l31 = 0.008

# At steady state (t -> +infinity), the derivatives are zero,
# which allows us to find the steady-state probabilities p_i.
# As derived in the explanation:
# p1 = p2 = p3
# p0 = (l10/l01) * p1
#
# Using the normalization condition p0 + p1 + p2 + p3 = 1, we get:
# p1 * (l10/l01 + 3) = 1
#
# The required sum is p0 + p1, which is derived as:
# p0 + p1 = (l10/l01 + 1) / (l10/l01 + 3)

# Calculate the ratio for convenience
ratio = l10 / l01

# Calculate the final result
result = (ratio + 1) / (ratio + 3)

# Output the final equation with all numbers and the result, as requested.
print(f"P_0(+inf) + P_1(+inf) = ({l10}/{l01} + 1) / ({l10}/{l01} + 3) = {result}")