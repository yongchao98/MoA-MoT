# The rarest noble gas on Earth is Radon (Rn) due to its radioactivity and short half-life.
# This script calculates its abundance as a percentage of Earth's total mass.

# Name of the rarest noble gas
gas_name = "Radon (Rn)"

# Estimated total mass of Radon on Earth in kilograms.
# This is derived from its average atmospheric abundance.
mass_gas = 3950.0  # kg

# Total mass of the Earth in kilograms.
mass_earth = 5.972e24  # kg

# Calculate the percentage of all terrestrial matter.
percentage = (mass_gas / mass_earth) * 100

print(f"The rarest noble gas on Earth is {gas_name}.")
print("To find its percentage of all terrestrial matter, we use the following equation:")
# We use .2e format specifier to represent the large number in a more readable scientific notation
print(f"({mass_gas} kg / {mass_earth:.3e} kg) * 100 = {percentage:.3e}%")
