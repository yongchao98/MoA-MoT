import math

# This script demonstrates the relationship between the radial wavevector k_r
# and the topological charge l for rotational propagation in Bessel-Gauss modes.
# The derived condition is k_r is proportional to sqrt(l).

# We can express this as k_r = C * sqrt(l), where C is a constant of proportionality.
# Let's assume a value for C for demonstration purposes.
C = 0.5  # Arbitrary proportionality constant

# A list of topological charges to test.
topological_charges = [1, 4, 9, 16]

print("Demonstrating the condition: k_r ∝ √l, or k_r = C * √l")
print(f"Using a proportionality constant C = {C}\n")

for l in topological_charges:
    # Calculate sqrt(l)
    sqrt_l = math.sqrt(l)
    # Calculate k_r based on the derived relationship
    k_r = C * sqrt_l
    
    # Output the relationship and the calculated values
    print(f"For topological charge l = {l}:")
    print(f"  The condition requires k_r = {C} * √{l}")
    print(f"  k_r = {C} * {sqrt_l:.2f}")
    print(f"  Resulting k_r = {k_r:.2f}\n")
