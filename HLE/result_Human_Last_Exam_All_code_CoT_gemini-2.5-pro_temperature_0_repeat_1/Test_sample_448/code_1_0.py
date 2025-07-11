import sys

# Define the constants needed for the calculation.
# The year 2002 does not change these fundamental physical constants.

# Mass of the Earth in kilograms.
mass_earth_kg = 5.972e24
# Mass of the Earth's atmosphere in kilograms.
mass_atmosphere_kg = 5.15e18
# Average molar mass of air in grams per mole.
molar_mass_air_g_mol = 28.97
# Molar mass of Radon's most stable isotope (Rn-222) in grams per mole.
molar_mass_radon_g_mol = 222
# Volume fraction of Radon in the atmosphere (unitless).
# Source: CRC Handbook of Chemistry and Physics (6 x 10^-18 percent -> 6e-20).
radon_volume_fraction = 6e-20

# Step 1: Calculate the mass of Radon in the atmosphere.
# The mass fraction is approximated by the volume fraction multiplied by the ratio of molar masses.
mass_radon_kg = mass_atmosphere_kg * (molar_mass_radon_g_mol / molar_mass_air_g_mol) * radon_volume_fraction

# Step 2: Calculate the percentage of Radon relative to the total mass of the Earth.
percentage_of_earth_mass = (mass_radon_kg / mass_earth_kg) * 100

# Step 3: Print the result, showing the numbers used in the final equation.
print("The rarest noble gas on Earth is Radon (Rn).")
print("Its percentage of all terrestrial matter is calculated below.")
print("-" * 50)
print(f"Mass of Radon (kg) = Mass of Atmosphere (kg) * (Molar Mass Radon / Molar Mass Air) * Volume Fraction")
print(f"Mass of Radon (kg) = {mass_atmosphere_kg:.3g} * ({molar_mass_radon_g_mol} / {molar_mass_air_g_mol}) * {radon_volume_fraction:.1g}")
print(f"Calculated Mass of Radon: {mass_radon_kg:.4g} kg")
print("-" * 50)
print(f"Percentage of Earth Mass = (Mass of Radon (kg) / Mass of Earth (kg)) * 100")
print(f"Percentage of Earth Mass = ({mass_radon_kg:.4g} / {mass_earth_kg:.4g}) * 100")
print(f"Final Result: {percentage_of_earth_mass:.4g}%")

# This is the final numerical answer for the tagging system.
# We need to redirect this to a different stream to avoid it being part of the main output.
# However, for this environment, we will just calculate it and use it for the final tag.
final_answer_value = f"{percentage_of_earth_mass:.4g}"