import sys
from io import StringIO

# The rarest noble gas on Earth is Radon (Rn).
# Its extreme rarity is due to it being radioactive with a short half-life.
# Its concentration in the atmosphere is highly variable but is estimated to be
# approximately 1 part in 10^21 (one sextillion). This is the best available
# figure to represent its percentage in terrestrial matter.

gas_name = "Radon"
# Representing 1 part per 10^21 as a fraction
concentration_as_fraction_numerator = 1
concentration_as_fraction_denominator = 1e21

# Convert the fraction to a percentage by multiplying by 100
percentage = (concentration_as_fraction_numerator / concentration_as_fraction_denominator) * 100

print(f"The rarest noble gas on Earth is {gas_name}.")
print(f"Its approximate concentration in the atmosphere is {concentration_as_fraction_numerator} part per {concentration_as_fraction_denominator:.0e}.")
print("To express this as a percentage, we perform the following calculation:")

# The final print statement shows the full equation as requested.
# Using scientific notation for clarity.
print(f"({concentration_as_fraction_numerator} / {concentration_as_fraction_denominator:.0e}) * 100 = {percentage:.1e}%")
