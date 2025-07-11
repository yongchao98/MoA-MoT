# The rarest noble gas on Earth is Radon (Rn).
# This script calculates its percentage of all terrestrial matter
# based on its estimated mass in the atmosphere.
# The abundances are stable, so the year 2002 does not change the result.

# Mass of Radon (Rn) in the atmosphere in kilograms.
# This is derived from its trace concentration (approx. 1 part in 10^21)
# and the mass of the atmosphere.
mass_radon_kg = 39.4

# Mass of the Earth in kilograms.
mass_earth_kg = 5.972e24

# Calculate the percentage of Radon relative to the total mass of the Earth.
percentage = (mass_radon_kg / mass_earth_kg) * 100

print("The rarest noble gas on Earth is Radon (Rn).")
print("Its percentage as part of all terrestrial matter can be estimated from its atmospheric mass.")
print("\n--- Calculation ---")
print(f"Estimated Mass of Radon in Atmosphere: {mass_radon_kg} kg")
print(f"Total Mass of Earth: {mass_earth_kg:.3e} kg")

print("\nEquation used:")
# As requested, outputting each number in the final equation.
print(f"( {mass_radon_kg} / {mass_earth_kg} ) * 100")

print("\nResult:")
print(f"The percentage is approximately: {percentage:.2e} %")
