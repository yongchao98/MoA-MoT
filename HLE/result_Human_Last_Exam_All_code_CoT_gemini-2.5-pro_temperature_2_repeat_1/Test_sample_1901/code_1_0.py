# Set a target transition temperature in Celsius, around room temperature.
target_temperature = 20.0

# Define the model parameters based on the design principles provided in the prompt.
# A pentyl chain (n=5) is the suggested starting point.
n_base = 5
# Assume a transition temperature for the base molecule (near room temperature).
T_base = 22.0  # Celsius
# Define the tuning factor (k): the change in transition temperature per carbon atom.
# The prompt states increasing chain length lowers the temperature, so k is positive.
# We will use a value of 4.0 C per carbon atom.
k = 4.0  # Celsius / carbon atom

# The model for the transition temperature T as a function of chain length n is:
# T(n) = T_base - k * (n - n_base)
# To find the 'n' for a target temperature, we rearrange the formula:
# n = n_base - (target_temperature - T_base) / k

# Calculate the ideal (potentially fractional) chain length for the target temperature.
ideal_n = n_base - (target_temperature - T_base) / k

# The number of carbons in the alkyl chain must be an integer, so we round to the nearest whole number.
calculated_n = int(round(ideal_n))

# Ensure the calculated chain length is a positive integer.
if calculated_n < 1:
    calculated_n = 1

# Calculate the final predicted transition temperature using the integer chain length.
predicted_T = T_base - k * (calculated_n - n_base)

# Determine the number of hydrogen atoms for the alkyl group (CnH2n+1).
num_hydrogens = 2 * calculated_n + 1

# Print the design recommendation and calculation steps.
print("Liquid Crystal Design Recommendation")
print("-" * 40)
print(f"Goal: Achieve a transition temperature near {target_temperature}°C.")
print(f"General Molecular Structure: CnH2n+1-Ph-CN")
print("-" * 40)
print("\nSuggested Molecule:")
print(f"Recommended alkyl chain length 'n' is: {calculated_n}")
print(f"The suggested molecule is 1-cyano-4-(C{calculated_n}H{num_hydrogens})benzene.")
print(f"Its predicted transition temperature is: {predicted_T:.1f}°C")

# Print the final equation with all numbers, as requested.
print("\nFinal Design Equation:")
print("T(n) = T_base - k * (n - n_base)")
print(f"{predicted_T:.1f} = {T_base} - {k} * ({calculated_n} - {n_base})")