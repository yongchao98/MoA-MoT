import math

# The number of terminal pairs given in the problem.
k = 10**6

# The best-known polynomial-time approximation algorithm for the Multicut problem
# on general graphs provides an approximation factor of O(log k).
# In the analysis of this algorithm, 'log' refers to the natural logarithm (ln).
# Let's calculate the value of log(k).
approximation_factor_val = math.log(k)

# Let's also calculate the value for the sqrt(log k) choice for comparison.
sqrt_log_k_val = math.sqrt(approximation_factor_val)

print(f"For k = {int(k)} terminal pairs:")
print(f"An O(log k) approximation gives a factor of ~{approximation_factor_val:.2f}.")
print(f"An O(sqrt(log k)) approximation would give a factor of ~{sqrt_log_k_val:.2f}.")

print("\nBased on established results, an O(log k) approximation is achievable, "
      "while a substantially better one (like O(sqrt(log k))) is not known to be possible in polynomial time.")
print("This means we can get an approximation α that satisfies the following relation:")

# The final output needs to show each number in the final equation.
# The relation is α ≤ log(k), which we can write with the calculated values.
print(f"α ≤ log({int(k)}) ≈ {approximation_factor_val:.1f}")