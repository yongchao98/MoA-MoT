# This script calculates the abundance of Xenon, the rarest stable noble gas,
# as a percentage of Earth's total mass. The values used are standard
# scientific constants that were the same in 2002.

# Step 1: Define the necessary constants.
# The mass of the Earth in kilograms.
mass_of_earth_kg = 5.972 * (10**24)

# The total mass of Xenon in Earth's atmosphere in kilograms.
# This is an excellent approximation for the total mass of Xenon on Earth.
mass_of_xenon_kg = 1.9 * (10**13)

# Step 2: Calculate the percentage.
# The formula is (mass of Xenon / mass of Earth) * 100.
percentage = (mass_of_xenon_kg / mass_of_earth_kg) * 100

# Step 3: Print the results and the components of the calculation.
print("The rarest stable noble gas on Earth is Xenon (Xe).")
print("Its abundance as a percentage of all terrestrial matter is calculated below.")
print("-" * 30)
print(f"Total mass of Xenon (kg): {mass_of_xenon_kg}")
print(f"Total mass of Earth (kg): {mass_of_earth_kg}")
print(f"Calculation: ({mass_of_xenon_kg} / {mass_of_earth_kg}) * 100")
print(f"Resulting Percentage: {percentage:.12f}%")
print("-" * 30)
print("Note: While Radon (Rn) is technically rarer, it is unstable and radioactive, so Xenon is considered the rarest stable noble gas.")