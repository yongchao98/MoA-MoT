# The rarest noble gas on Earth is Radon (Rn).
# This script calculates its estimated abundance as a percentage of all terrestrial matter.

# Step 1: Define the known and estimated values.
# The year 2002 does not change these fundamental physical quantities.
gas_name = "Radon (Rn)"
# Estimated total mass of Radon in the Earth's crust in kilograms.
# This is derived from the equilibrium concentration with its parent element, Uranium.
mass_radon_kg = 1.3e5
# Total mass of the Earth in kilograms.
mass_earth_kg = 5.972e24

# Step 2: Calculate the percentage.
# Percentage = (part / whole) * 100
percentage = (mass_radon_kg / mass_earth_kg) * 100

# Step 3: Print the result, showing the numbers used in the calculation.
print(f"The rarest noble gas on Earth is {gas_name}.")
print(f"Its estimated total mass on Earth is approximately {mass_radon_kg:.1e} kg.")
print(f"The total mass of the Earth is approximately {mass_earth_kg:.3e} kg.")
print(f"The percentage abundance is calculated as: ({mass_radon_kg} / {mass_earth_kg}) * 100")
print(f"The percentage of Radon in all terrestrial matter is approximately: {percentage:.20f} %")

# Final answer variable for direct output
final_answer_value = percentage