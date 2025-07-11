import math

def calculate_force_isolated_polymer(x, n, l, E0):
  """
  Calculates the restoring force for a thermally isolated freely jointed chain.
  This is the force exerted by the polymer on its ends.

  Args:
    x: The end-to-end separation of the polymer.
    n: The number of segments (mass points).
    l: The length of each segment.
    E0: The kinetic energy of the polymer at zero extension, E(0).

  Returns:
    The restoring force F(x). The force is negative, indicating it's attractive.
  """
  # Check if the Gaussian approximation is likely valid.
  if abs(x) > n * l:
    print("Warning: The extension x is large; this formula may not be accurate.")

  # The derived force law is F(x) = - (2 * E(0) * x / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))
  n_sq_l_sq = (n**2) * (l**2)
  if n_sq_l_sq == 0:
      return 0

  numerator = 2 * E0 * x
  exponent_term = (x**2) / n_sq_l_sq
  
  force = -(numerator / n_sq_l_sq) * math.exp(exponent_term)
  return force

# The primary task is to find the force law. The law is a symbolic formula.
# We print the derived formula with its components, as requested.

print("The derived force law for a thermally isolated polymer is given by the formula:")
print("F(x) = - ( (2) * (E(0)) * (x) / ( (n)^2 * (l)^2 ) ) * exp( (x)^2 / ( (n)^2 * (l)^2 ) )")
print("\nWhere:")
print("  F(x): Force as a function of extension x")
print("  E(0): Kinetic energy at zero extension (x=0)")
print("  x:    Separation of the polymer ends")
print("  n:    Number of mass points (segments)")
print("  l:    Length of a segment")

# --- Example Calculation ---
# You can uncomment the following lines to run an example calculation.
# print("\n--- Example Calculation ---")
# x_val = 1.0e-8   # meters
# n_val = 1000     # segments
# l_val = 1.0e-9   # meters per segment
# E0_val = 2.07e-21 # Joules (corresponds to T~100K for n=1000)
# F_calculated = calculate_force_isolated_polymer(x_val, n_val, l_val, E0_val)
# print(f"For E(0) = {E0_val:.2e} J, n = {n_val}, l = {l_val:.2e} m, and x = {x_val:.2e} m,")
# print(f"the calculated force is F = {F_calculated:.4e} N")
