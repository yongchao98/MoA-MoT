import math

# The number of terminal pairs in the Multicut problem instance.
k = 10**6

# In the field of approximation algorithms, it is a well-established result
# that the Multicut problem on general graphs has a polynomial-time
# O(log k)-approximation algorithm. Furthermore, this is considered a tight bound,
# as it's known to be NP-hard to approximate the problem to a factor better than c*log(k)
# for some constant c.

# We will now calculate the values for the approximation ratios mentioned
# in the answer choices for k = 10^6.

# Approximation ratio from choice B: sqrt(log k)
sqrt_log_k = math.sqrt(math.log(k))

# Approximation ratio from choice C: log k
log_k = math.log(k)

# Print the analysis.
print("Analysis for Multicut with k = 1,000,000 terminal pairs:")
print(f"An O(sqrt(log k)) approximation would be: sqrt(log({k:.0f})) ~= {sqrt_log_k:.1f}")
print(f"An O(log k) approximation would be: log({k:.0f}) ~= {log_k:.1f}")
print("\nAccording to known results, an O(log k) approximation is achievable in polynomial time.")
print("The hardness results show that we cannot do significantly better, so an O(sqrt(log k)) approximation is not known to be possible.")
print("This means the best known approximation factor is alpha <= log(k).")
print("\nFinal equation for the approximation factor alpha:")
# The final result prints each number in the equation.
print(f"alpha <= log({k}) ~= {log_k:.1f}")
