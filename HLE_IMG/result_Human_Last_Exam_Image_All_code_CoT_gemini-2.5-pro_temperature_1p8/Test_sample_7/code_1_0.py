# Define the known voltage values for the two plateaus.
V1 = 0.09  # Voltage of the first plateau (50%-100% SOC)
V2 = 0.13  # Voltage of the second plateau (20%-50% SOC)

# Calculate the voltage difference, which corresponds to the chemical potential difference term.
voltage_difference = V2 - V1

# The problem asks for a formula for the second plateau (V2).
# We derived the relationship V2 = V1 + (μ_1 - μ_2) / e.
# Let's print this formula and substitute the known numerical values.

print("The relationship between the two voltage plateaus is given by:")
print("V_2 - V_1 = (μ_1 - μ_2) / e")
print("\nRearranging to find a formula for the second plateau (V_2):")
print("V_2 = V_1 + (μ_1 - μ_2) / e")
print("\nSubstituting the numerical values given in the problem:")
# The prompt requires outputting each number in the final equation.
print(f"{V2:.2f} V = {V1:.2f} V + ({voltage_difference:.2f} V)")
print("\nwhere the term (μ_1 - μ_2) / e is equal to the voltage difference between the plateaus.")
