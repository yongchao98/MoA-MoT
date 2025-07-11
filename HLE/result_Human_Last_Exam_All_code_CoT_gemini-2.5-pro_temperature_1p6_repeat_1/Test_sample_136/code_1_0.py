import math

# Step 1: Define the probability u based on the theoretical derivation.
# u is the probability that a walk starting from the origin 0 hits x_0
# before returning to 0.
# u = 1/2 * (h(c) + h(u)), where c is a common neighbor and u is an uncommon one.
# h(c) = P_c(hit x_0 before 0) = 1/2
# h(u) = P_u(hit x_0 before 0) = 1/2 - 1/pi

pi = math.pi
h_c = 1/2
h_u = 1/2 - 1/pi

u = 1/2 * (h_c + h_u)

# Step 2: The limit L is given by the formula L = 1 - u^2.
L = 1 - u**2

# Step 3: Print the equation with the intermediate values.
print(f"The probability u is calculated as:")
print(f"u = 1/2 * (P_common_neighbor(hit x_0 first) + P_uncommon_neighbor(hit x_0 first))")
print(f"u = 1/2 * ({h_c} + ({h_u:.4f}))")
print(f"u = 1/2 * (({h_c + h_u:.4f}))")
print(f"u = {u:.4f}")
print("\nThe final limit L is calculated as:")
print(f"L = 1 - u^2")
print(f"L = 1 - ({u:.4f})^2")
print(f"L = 1 - {u**2:.4f}")
print(f"L = {L:.4f}")

# The symbolic result is 3/4 + 1/(2*pi) - 1/(4*pi^2)
symbolic_L = 3/4 + 1/(2*pi) - 1/(4*pi**2)
print("\nFor verification, the symbolic calculation gives:")
print(f"L = 3/4 + 1/(2*pi) - 1/(4*pi^2) = {symbolic_L:.4f}")

print(f"\nFinal numerical answer:")
print(L)
