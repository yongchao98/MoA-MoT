import math

# The problem asks for the best polynomial-time approximation factor for the Multicut Problem.
# The number of terminal pairs is k.
k = 10**6

# The seminal work by Garg, Vazirani, and Yannakakis provides an algorithm with an
# approximation guarantee of O(log k). This is the best-known result for general graphs.
# In this context, "log" refers to the natural logarithm (ln).

# We will now calculate the value of log(k) for the given k.
approximation_factor = math.log(k)

# The final answer choice C states that the approximation is around log(k), which is approximately 13.8.
# Our code will verify this calculation.
# As requested, we will print the final equation with each number.
print("The best-known polynomial-time approximation factor alpha for the Multicut problem is O(log k).")
print("For k = 1,000,000, the approximation factor is calculated as:")
print(f"alpha <= log({int(k)})")
print(f"alpha <= {approximation_factor:.1f}")
