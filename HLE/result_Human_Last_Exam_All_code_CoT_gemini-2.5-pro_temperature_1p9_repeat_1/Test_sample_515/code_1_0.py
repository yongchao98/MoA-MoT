import math

# --- Constants ---
# Value of pi from the math library
pi = math.pi
# Euler-Mascheroni constant
gamma = 0.5772156649

# --- Start point potential ---
# The starting point is x0 = (0,1). Its potential is a(0,1).
# By symmetry, a(0,1) = a(1,0). The exact value is 4/pi - 1.
a_x0 = 4/pi - 1

# --- Target point potential ---
# The target set A is the neighbors of z = (3600,0).
# Since the distance is large, we use the asymptotic formula for a(n,0):
# a(n,0) ~= (2/pi) * (ln(n) + gamma + (3/2)*ln(2))
# We approximate the potential of the set A by the potential at its center z.
n = 3600
log_n = math.log(n)
log_2 = math.log(2)
a_z_approx = (2/pi) * (log_n + gamma + 1.5 * log_2)

# --- Probability calculation ---
# The probability is approximated by the ratio of potentials a(x0)/a(z)
probability = a_x0 / a_z_approx

# --- Output the result ---
# Print the components of the calculation for clarity.
print("Calculation of the probability P that the conditioned random walk hits the target set.")
print("The probability is approximated by the ratio of potentials P ≈ a(x0) / a(z).")
print(f"Starting point x0 = (0,1)")
print(f"Potential at starting point a(x0) = 4/π - 1 = {a_x0}")
print(f"Target point z = (3600,0)")
print(f"Asymptotic potential at target point a(z) ≈ (2/π) * (ln({n}) + γ + 1.5*ln(2)) = {a_z_approx}")
print("\nFinal Probability Calculation:")
# Display the full equation being computed.
# We represent the equation with numbers rounded for display, but use full precision for the calculation.
print(f"P ≈ {a_x0:.4f} / {a_z_approx:.4f}")
print(f"P ≈ {probability:.5f}")

# Give the answer with two significant digits as requested
print(f"\nThe approximate answer with two significant digits is: {probability:.2g}")