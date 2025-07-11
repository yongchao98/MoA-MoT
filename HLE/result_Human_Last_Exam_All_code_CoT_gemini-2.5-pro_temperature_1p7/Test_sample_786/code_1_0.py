import math

# The number of terminal pairs in the Multicut problem.
k = 10**6

# The best-known and asymptotically tight polynomial-time approximation factor for the Multicut
# problem is O(log k). For the specific value, we use the natural logarithm (ln) as
# indicated by the problem's options.
# We will calculate the value for this approximation guarantee.
approximation_factor_val = math.log(k)

# We present the result in the format of the final equation as requested.
print(f"The best polynomial-time approximation factor α for the Multicut problem with k terminal pairs is bounded by O(log k).")
print(f"For k = {int(k)}, the equation for the approximation factor is:")
print(f"α ≤ log({int(k)})")
print(f"Calculating the value:")
print(f"α ≤ {approximation_factor_val:.2f}")
print("\nThis matches option C.")
