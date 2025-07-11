import sys

# Using float('inf') for Radon as its concentration is exceptionally low and variable,
# often cited as 'trace amounts' or in atoms/m^3 rather than a stable percentage.
# For practical comparison, we can use an extremely small number to represent its percentage,
# such as 1e-19 %, which ensures it will be the minimum.
# Data is based on atmospheric composition by volume.
noble_gases = {
    "Helium": 0.000524,
    "Neon": 0.001818,
    "Argon": 0.934,
    "Krypton": 0.000114,
    "Xenon": 0.0000087,
    "Radon": 1e-19  # ~1 part per 1,000,000,000,000,000,000,000 or 1e-19 %
}

print("Noble Gas Abundance on Earth (as % of atmosphere):")
for gas, percentage in noble_gases.items():
    # Use scientific notation for very small numbers for better readability
    if percentage < 1e-6:
        print(f"{gas} < {percentage:.2e}%")
    else:
        print(f"{gas} < {percentage}%")

# Find the rarest gas by finding the minimum percentage in the dictionary
rarest_gas_name = min(noble_gases, key=noble_gases.get)
rarest_gas_percentage = noble_gases[rarest_gas_name]

print("\n-----------------------------------------------------")
print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
sys.stdout.write("Its abundance as a percentage of all terrestrial matter is < ")
# sys.stdout.write is used to avoid a newline to construct the final output line
# The final line will show the rarest gas's percentage.
print(f"{rarest_gas_percentage:.2e}%")
print("-----------------------------------------------------")