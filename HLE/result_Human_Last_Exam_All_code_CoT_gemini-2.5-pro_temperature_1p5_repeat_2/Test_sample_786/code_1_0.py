import math

# The number of terminal pairs in the Multicut problem.
k = 10**6

# The best-known polynomial-time approximation ratio for Multicut is O(log k).
# In theoretical computer science, 'log' typically refers to the natural logarithm (ln)
# unless specified otherwise. We calculate this value for the given k.
approximation_ratio = math.log(k)

print("For the Multicut Problem:")
print(f"The number of terminal pairs is k = {k}.")
print("The best achievable polynomial-time approximation ratio is O(log k).")
print("\nCalculating the value for k = 10^6:")
# The problem statement requires printing all numbers in the final equation.
print(f"alpha <= log({int(k)}) approx {approximation_ratio:.1f}")
print("\nThis result aligns with option C.")
