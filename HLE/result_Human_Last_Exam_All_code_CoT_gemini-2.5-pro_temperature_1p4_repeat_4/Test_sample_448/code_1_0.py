# The rarest noble gas on Earth is Radon (Rn).
# Its extreme rarity is due to its radioactivity and short half-life (the most stable isotope, Radon-222, has a half-life of about 3.8 days).
# Its concentration in the atmosphere is minuscule, on the order of 1 part per 10^21.

# Define the gas and its abundance details.
gas_name = "Radon"

# Abundance as a fraction (1 part in 10^21)
fractional_abundance = 1e-21

# Convert the fractional abundance to a percentage.
# Percentage = Fraction * 100
percentage_abundance = fractional_abundance * 100

# Print the final result, showing the numbers used.
print(f"The rarest noble gas on Earth as a percentage of atmospheric matter is {gas_name}.")
print(f"Its abundance is approximately 1 part per 10^21.")
print(f"To express this as a percentage, we calculate (1 / 1e21) * 100.")
print(f"The resulting percentage is: {percentage_abundance:.1e}%")
