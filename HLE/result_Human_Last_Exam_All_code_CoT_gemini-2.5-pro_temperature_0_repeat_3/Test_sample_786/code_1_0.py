import math

# Number of terminal pairs
k = 10**6

# Calculate log k (natural logarithm)
log_k = math.log(k)

# Calculate sqrt(log k)
sqrt_log_k = math.sqrt(log_k)

# The best known approximation ratio is O(log k).
# Let's print the values relevant to the options.
print(f"The number of terminal pairs is k = {int(k)}")
print(f"The approximation factor based on log k is approximately: log({int(k)}) = {log_k:.1f}")
print(f"The value for the incorrect sqrt(log k) factor would be: sqrt(log({int(k)})) = {sqrt_log_k:.1f}")
