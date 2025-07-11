import math

# Step 1: The problem provides the slice rank formula SR = (3/2^K)^n * e^(o(n)).
# Step 2: We derived a lower bound for the slice rank, SR >= 2^n + 1.
# Step 3: We assume the asymptotic behavior of the slice rank is determined by this lower bound.
# This leads to the equation: 3 / (2^K) = 2.
# We solve for K.

# 3 / (2**K) = 2  => 2**K = 3/2
# K = log2(3/2)

K = math.log2(3) - math.log2(2)

print("The final equation is 3 / (2^K) = 2.")
print("Solving for K, we get 2^K = 3/2.")
print("This gives K = log2(3/2) = log2(3) - 1.")
print(f"The value of K is: {K}")
