import math

# Number of terminal pairs
k = 10**6

# Calculate log k (natural logarithm)
log_k = math.log(k)

# Calculate sqrt(log k)
sqrt_log_k = math.sqrt(log_k)

print(f"For k = {int(k)}:")
print(f"The value for a 'log k' approximation is: log({int(k)}) = {log_k:.1f}")
print(f"The value for a 'sqrt(log k)' approximation is: sqrt(log({int(k)})) = {sqrt_log_k:.1f}")
