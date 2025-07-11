import math

# Number of terminal pairs
k = 10**6

# Calculate log(k)
log_k = math.log(k)

# Calculate sqrt(log k)
sqrt_log_k = math.sqrt(log_k)

# The problem asks for an approximation of alpha <= log(k)
# We will print the equation and the result.
# Remember in the final code you still need to output each number in the final equation!
print(f"The number of terminal pairs k is {k}.")
print(f"The approximation factor based on log(k) is calculated as:")
print(f"log({k}) = {log_k:.1f}")
print(f"An alternative but unachievable approximation factor would be:")
print(f"sqrt(log({k})) = {sqrt_log_k:.1f}")