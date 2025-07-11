# The problem is to find the value of sum_{k=1 to 6} k * n_k based on scaling laws
# for magnetic field noise from a conducting slab.

# Exponents for the zero-frequency scaling:
# S_B(omega -> 0) proportional to sigma^(n_1) * T^(n_2) * z^(n_3)
n1 = 2  # Proportional to sigma^2
n2 = 1  # Proportional to T^1
n3 = -2 # Proportional to z^(-2)

# Exponents for the frequency-dependent scaling S_B(omega) proportional to omega^(n_k)
# in different regimes.
n4 = 2    # For omega << 1/(sigma*z*t), S_B scales as omega^2
n5 = 0    # For 1/(sigma*z*t) << omega << 1/(sigma*t^2), S_B is constant
n6 = -0.5 # For omega >> 1/(sigma*t^2), S_B scales as omega^(-1/2)

# Store the exponents in a list for easier access
n = [n1, n2, n3, n4, n5, n6]

# Calculate the sum
total_sum = 0
for k in range(1, 7):
    # The list is 0-indexed, so we use k-1
    term = k * n[k-1]
    total_sum += term
    print(f"k={k}, n_{k}={n[k-1]}, k*n_{k} = {term}")

print("\nThe final equation is:")
print(f"1*({n1}) + 2*({n2}) + 3*({n3}) + 4*({n4}) + 5*({n5}) + 6*({n6}) = {total_sum}")

print(f"\nThe final value of the sum is: {total_sum}")