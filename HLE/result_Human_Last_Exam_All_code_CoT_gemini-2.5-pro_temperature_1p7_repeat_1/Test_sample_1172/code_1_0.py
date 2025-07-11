# This script will print the derived formula for the change in mutual inductance.

# Define the terms of the formula as strings for clear output.
# This represents the approximate mutual inductance of the bare circuits (M1)
term_M1 = "-(mu_0 * h**2) / (2 * pi * d**2)"

# This represents the scaling factor introduced by the concentrator
term_scaling = "((R2/R1)**2 - 1)"

print("The expression for the change in mutual inductance per unit length (Delta_M),")
print("in the limit where the circuit separation d is much greater than the wire spacing h, is:")
print("\nDelta_M = M1 * [ (R2/R1)**2 - 1 ]")
print("\nSubstituting the expression for M1, we get:")
# The final equation is constructed by combining the different terms.
print(f"\nDelta_M = ({term_M1}) * ({term_scaling})")

print("\nWhere:")
print("  mu_0: Permeability of free space")
print("  h:    Separation distance between wires in each circuit")
print("  d:    Separation distance between the two circuits")
print("  R1:   Inner radius of the concentrator shell")
print("  R2:   Outer radius of the concentrator shell")

final_expression = f"- (mu_0 * h**2 / (2 * pi * d**2)) * ((R2/R1)**2 - 1)"
print(f"\nFinal formatted expression:\n\n{final_expression}")