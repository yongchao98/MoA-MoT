import math

# The problem is to find K in the expression for the slice rank of a tensor.
# The slice rank is given as (3 / 2**K)**n * exp(o(n)).
# This means the asymptotic slice rank is gamma = 3 / 2**K.

# We determined that the asymptotic slice rank of the tensor is 2.
# So, we need to solve the equation: gamma = 2.
gamma = 2
val_3 = 3
val_2 = 2

# The equation is 3 / (2**K) = 2
# 3 = 2 * (2**K)
# 3 / 2 = 2**K
# K = log2(3/2)
K = math.log2(val_3 / val_2)

print("The final equation to solve for K is:")
print(f"{val_3} / ({val_2}**K) = {val_2}")
print(f"Solving for K gives K = log2(3/2) = log2(3) - 1")
print(f"The value of K is: {K}")