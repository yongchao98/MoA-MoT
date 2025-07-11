import math

# The number of terminal pairs
k = 10**6

# The best-known approximation factor for the Multicut problem is O(log k).
# The answer choice C suggests using the natural logarithm (ln).
# We calculate the value of log k for k = 10^6.
approximation_factor_val = math.log(k)

print("The best known polynomial-time approximation for the Multicut problem is alpha <= O(log k).")
print(f"For k = {k}, we calculate the value of the approximation factor:")
print(f"alpha <= log({k}) approx {approximation_factor_val:.5f}")
