import math

# The number of terminal pairs
k = 10**6

# For the general Multicut problem, the best-known polynomial-time
# approximation algorithm (by Garg, Vazirani, and Yannakakis)
# achieves an approximation factor of O(log k).

# We calculate the value of log k. The problem's options suggest
# using the natural logarithm (base e).
log_k = math.log(k)

print(f"The Multicut problem has k = {int(k)} terminal pairs.")
print("The best-known polynomial-time approximation algorithm provides an approximation ratio (alpha) related to log(k).")
print("\nThe approximation relationship is: alpha <= log(k)")
print(f"Substituting the value of k = {int(k)}:")
print(f"alpha <= log({int(k)})")
print(f"alpha <= {log_k:.1f}")
