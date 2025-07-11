import sys

# This script performs a unit conversion for the patient's weight and states the oxygen saturation.
# The values are taken from the problem description.
pounds = 12
ounces = 1
oxygen_saturation = 89

# Constants for conversion
OUNCES_IN_A_POUND = 16
POUNDS_TO_KG = 0.453592

# Perform the calculation
total_pounds = pounds + (ounces / OUNCES_IN_A_POUND)
total_kg = total_pounds * POUNDS_TO_KG

# Per the user's request, print out each number in the final equation.
# The equation demonstrates the conversion of the baby's weight from pounds and ounces to kilograms.
# We will format the output to be clear and show all the original numbers.
print(f"Patient weight: {pounds} lb {ounces} oz")
print(f"Oxygen Saturation: {oxygen_saturation}%")
print("\n--- Weight Conversion Equation ---")
# The following line explicitly shows each number used in the final calculation
print(f"Calculation: ({pounds} + {ounces}/{OUNCES_IN_A_POUND}) * {POUNDS_TO_KG} = {total_kg:.2f} kg")