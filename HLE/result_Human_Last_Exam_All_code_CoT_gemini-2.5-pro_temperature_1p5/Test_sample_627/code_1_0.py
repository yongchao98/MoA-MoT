import numpy as np

# This script calculates the upper bound for the braid index of the 
# three-twist knot (6_1) using Vogel's algorithm.

# Step 1: Define the Alexander polynomial for the 6_1 knot.
# The standard (symmetrized) Alexander polynomial for 6_1 is 2t^2 - 3t + 2.
# We can represent this as a list of coefficients [2, -3, 2] for powers t^0, t^1, t^2.
# Using numpy's polynomial class for easier calculations, coefficients are [c0, c1, c2, ...].
delta_k_coeffs = [2, -3, 2]
delta_k = np.poly1d(delta_k_coeffs[::-1]) # numpy wants coefficients from highest power

print("Step-by-step calculation for the braid index upper bound of the three-twist knot (6_1):")
print("-" * 70)

print(f"1. The Alexander polynomial for the 6_1 knot is: Δ_K(t) = {delta_k}\n")

# Step 2: For an alternating knot, Vogel's polynomial P_K(t) is given by:
# P_K(t) = (t+1) * Δ_K(-t) - (t-1)
# The braid index bound is b(K) <= 1 + deg(P_K(t)).
print("2. For an alternating knot, Vogel's polynomial P_K(t) is calculated as:")
print("   P_K(t) = (t+1) * Δ_K(-t) - (t-1)\n")

# Step 3: Calculate Δ_K(-t)
# This is equivalent to multiplying coefficients of odd powers by -1.
delta_k_minus_t = delta_k.copy()
# In polynomial p(x) = c0 + c1*x + c2*x^2 + ..., p(-x) = c0 - c1*x + c2*x^2 - ...
# The coefficients of numpy polynomial are [c_n, c_{n-1}, ..., c0]
# So we negate coeffs at indices corresponding to odd powers: n-1, n-3, ...
for i in range(len(delta_k_minus_t.coeffs)):
    power = len(delta_k_minus_t.coeffs) - 1 - i
    if power % 2 != 0:
        delta_k_minus_t.coeffs[i] *= -1

print(f"3. First, we find Δ_K(-t). Substituting t with -t in Δ_K(t):")
print(f"   Δ_K(-t) = {delta_k_minus_t}\n")

# Step 4: Calculate (t+1) * Δ_K(-t)
t_plus_1 = np.poly1d([1, 1])
term1 = t_plus_1 * delta_k_minus_t

print("4. Next, multiply Δ_K(-t) by (t+1):")
print(f"   (t+1) * ({delta_k_minus_t}) = {term1}\n")

# Step 5: Calculate P_K(t) by subtracting (t-1)
t_minus_1 = np.poly1d([1, -1])
p_k = term1 - t_minus_1

print("5. Now, we subtract (t-1) to get the full Vogel polynomial P_K(t):")
print(f"   P_K(t) = ({term1}) - ({t_minus_1})")
print(f"   P_K(t) = {p_k}\n")

# Step 6: Find the degree of P_K(t)
degree = p_k.order
print(f"6. The degree of P_K(t) is the highest power of t, which is: {degree}\n")

# Step 7: Calculate the upper bound for the braid index
upper_bound = 1 + degree
print("7. Finally, the upper bound for the braid index b(K) is 1 + deg(P_K(t)):")
print(f"   b(K) <= 1 + {degree}")
print(f"   b(K) <= {upper_bound}\n")

print(f"The upper bound for the braid index of the three-twist knot is {upper_bound}.")
