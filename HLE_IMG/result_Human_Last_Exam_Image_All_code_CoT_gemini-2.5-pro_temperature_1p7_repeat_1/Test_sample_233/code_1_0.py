# The runestone ID consists of a province code and a number.
# The province for this runestone is Uppland, represented by 'U'.
# The number for this specific stone is 778.

# As per the instructions, let's represent the number part using an equation.
province_code = "U"
hundreds = 700
tens = 70
units = 8

# Calculate the final number from its components.
id_number = hundreds + tens + units

# Print the final ID, showing the equation as requested.
print("The ID of the Ingvar runestone is calculated as follows:")
print(f"Province Code: '{province_code}'")
print(f"Number Calculation: {hundreds} + {tens} + {units} = {id_number}")
print("---")
print(f"Final Runestone ID: {province_code}{id_number}")