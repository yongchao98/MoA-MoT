import sys

# Step 1: Define the exponents based on the physical analysis.
# n1, n2, n3 from the zero-frequency scaling: S_B ~ T * sigma * t * z^(-4)
n1 = 1
n2 = 1
n3 = -4

# n4, n5, n6 from the frequency spectrum scaling
# Regime omega << 1/(sigma*z*t) -> transparent, non-shielding -> white noise
n4 = 0
# Regime 1/(sigma*z*t) << omega << 1/(sigma*t^2) -> transparent, shielding -> S_B ~ omega^(-2)
n5 = -2
# Regime omega >> 1/(sigma*t^2) -> opaque, shielding -> S_B ~ omega^(-3/2)
n6 = -1.5

exponents = [n1, n2, n3, n4, n5, n6]
total_sum = 0
sum_expression_parts = []

# Step 2: Calculate the sum and create the expression string
for k in range(1, 7):
    n_k = exponents[k - 1]
    term = k * n_k
    total_sum += term
    sum_expression_parts.append(f"({k})*({n_k})")

sum_expression = " + ".join(sum_expression_parts)

# Step 3: Print the full equation and the final answer.
# Using sys.stdout.write to ensure no extra characters are printed,
# as per the implied instruction to not alter the final output format.
output_string = f"The sum is calculated as:\n{sum_expression} = {total_sum}\n"
sys.stdout.write(output_string)

# The final result in the requested format
sys.stdout.write(f"\n<<<{total_sum}>>>\n")
