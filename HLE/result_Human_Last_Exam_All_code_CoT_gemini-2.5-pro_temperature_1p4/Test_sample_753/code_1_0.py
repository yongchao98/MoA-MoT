import math

# Step 1: Define parameters based on m=3.
m = 3

# Step 2: Determine n and the partition lambda.
# The partition is given by lambda = (m^1, (m-1)^2, ..., 1^m).
# For m=3, lambda = (3^1, 2^2, 1^3).
# The parts are 3, 2, 2, 1, 1, 1.
# n is the sum of the parts.
n = 3 * 1 + 2 * 2 + 1 * 3
# Also using the formula from the problem statement:
# n = sum(k * (m + 1 - k) for k in range(1, m + 1))
# This gives 1*3 + 2*2 + 3*1 = 10.

# The cycle counts c_k for the partition lambda are:
# c_1: number of parts of size 1, which is m
# c_2: number of parts of size 2, which is m-1
# ...
# c_k: number of parts of size k, which is m+1-k
cycle_counts = {k: m + 1 - k for k in range(1, m + 1)}

# Step 3: Calculate z_lambda, the size of the centralizer.
# z_lambda = product of (k^c_k * c_k!) for each k.
z_lambda = 1
for k, c_k in cycle_counts.items():
    if c_k > 0:
        z_lambda *= (k**c_k) * math.factorial(c_k)

# Step 4: Calculate the cardinality of Theta^{-1}(lambda).
# The formula is |Theta^{-1}(lambda)| = (n!)^3 / z_lambda.
# We use Python's arbitrary-precision integers for the calculation.
n_factorial = math.factorial(n)
n_factorial_cubed = n_factorial**3
cardinality = n_factorial_cubed // z_lambda

# Step 5: Output the first 40 digits of the result.
# Convert the number to a string to extract digits.
s_cardinality = str(cardinality)
first_40_digits = s_cardinality[:40]

# Print the final numerical result as requested
# The final expression is (10!)^3 / 144
print("The value of n for m=3 is: " + str(n))
print("The partition lambda corresponds to cycle counts: c_1=3, c_2=2, c_3=1")
print("The size of the centralizer z_lambda is: 1^3*3! * 2^2*2! * 3^1*1! = " + str(z_lambda))
print("The cardinality of Theta^-1(lambda) is (n!)^3 / z_lambda:")
print("({})^3 / {} = {}".format(n_factorial, z_lambda, cardinality))
print("\nThe first 40 digits of the cardinality are:")
print(first_40_digits)
