# The rarest noble gas is Radon (Rn) due to its radioactive nature and short half-life.
# This script calculates its approximate percentage of all terrestrial matter.

# Step 1: Define the estimated total mass of Radon on Earth in kilograms.
# This value is derived from the steady-state equilibrium of Radon-222
# with its parent elements in the Earth's crust and mantle.
total_radon_mass_kg = 2.2e5

# Step 2: Define the total mass of the Earth in kilograms.
total_earth_mass_kg = 5.972e24

# Step 3: Calculate the abundance as a percentage.
# The formula is (Mass of Radon / Mass of Earth) * 100.
percentage = (total_radon_mass_kg / total_earth_mass_kg) * 100

# Step 4: Print the result and the equation used.
print("The rarest noble gas on Earth is Radon (Rn).")
print(f"Its estimated total mass on Earth is approximately {total_radon_mass_kg:e} kg.")
print(f"The total mass of the Earth is approximately {total_earth_mass_kg:e} kg.")
print("\nThe calculation for its percentage of all terrestrial matter is:")
print(f"({total_radon_mass_kg} / {total_earth_mass_kg}) * 100 = {percentage:e}%")
