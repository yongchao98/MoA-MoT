import collections

# Step 1 & 2: Define a dictionary with the noble gases and their atmospheric abundance in percent.
# The abundance of Radon is extremely low and variable; this is an estimated average value.
noble_gas_abundances = collections.OrderedDict({
    "Helium": 0.000524,
    "Neon": 0.001818,
    "Argon": 0.934,
    "Krypton": 0.000114,
    "Xenon": 0.0000087,
    "Radon": 6e-18,
})

# Step 3: Find the rarest gas (the one with the minimum abundance).
rarest_gas_name = min(noble_gas_abundances, key=noble_gas_abundances.get)
rarest_gas_percentage = noble_gas_abundances[rarest_gas_name]

# Step 4: Format the search for the minimum value as an equation,
# fulfilling the "output each number in the final equation" requirement.
all_percentages = list(noble_gas_abundances.values())
equation_str = f"min({', '.join(map(str, all_percentages))}) = {rarest_gas_percentage}"

print("To find the rarest noble gas, we compare their atmospheric percentages:")
print(equation_str)
print("\n---\n")

# Step 5: Print the conclusion.
print(f"The minimum percentage in the comparison is {rarest_gas_percentage}%.")
print(f"This percentage corresponds to the noble gas: {rarest_gas_name}.")
print(f"\nTherefore, the rarest noble gas on Earth as a percentage is {rarest_gas_name} with approximately {rarest_gas_percentage}% abundance in the atmosphere.")
