import math

# Step 1: The problem is to find K from the slice rank expression.
# The slice rank is given as (3/2^K)^n * e^(o(n)).
# Our analysis shows that the slice rank of the tensor is exactly 2^n.

# Step 2: Equate the two expressions for the slice rank.
# 2^n = (3 / (2^K))^n * e^(o(n))
# Taking the n-th root of both sides and letting n -> infinity, we get:
# 2 = 3 / (2^K)

# Step 3: Solve the equation for K.
# 2 * 2^K = 3
# 2^(K+1) = 3
# K + 1 = log2(3)
# K = log2(3) - 1

# The code below calculates and prints the value of K.

# The numbers in the final equation K = log2(3) - 1 are K, 2, 3, and 1.
log_3_base_2 = math.log2(3)
one = 1
K = log_3_base_2 - one

print("The final equation for K is: K = log\u2082(3) - 1")
print(f"We compute the values: log\u2082(3) = {log_3_base_2}, and 1 = {one}")
print(f"So, K = {log_3_base_2} - {one}")
print(f"The value of K is: {K}")
