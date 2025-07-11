import math

# Number of terminal pairs
k = 10**6

# Calculate log(k) using the natural logarithm
log_k = math.log(k)

# Calculate sqrt(log(k))
sqrt_log_k = math.sqrt(log_k)

print("The problem considers k terminal pairs.")
print(f"k = {int(k)}")
print("\nThe best-known approximation factor is O(log k).")
print(f"log(k) = log({int(k)}) ≈ {log_k:.1f}")
print("\nThis rules out better approximations like O(sqrt(log k)).")
print(f"sqrt(log(k)) = sqrt({log_k:.1f}) ≈ {sqrt_log_k:.1f}")
print("\nTherefore, we can get an alpha <= log k approximation, but not an alpha <= sqrt(log k) approximation.")
