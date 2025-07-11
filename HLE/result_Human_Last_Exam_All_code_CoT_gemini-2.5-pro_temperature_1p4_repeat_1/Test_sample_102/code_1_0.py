import math

# Define the given equivalent values
cr_eq = 39
ni_eq = 29

# Constants from the empirical formula for ferrite level calculation
cr_factor = 3.34
ni_factor = 2.46
constant_term = 28.6

# Calculate the ferrite level using the formula:
# Ferrite Level = (cr_factor * Cr_eq) - (ni_factor * Ni_eq) - constant_term
ferrite_level = (cr_factor * cr_eq) - (ni_factor * ni_eq) - constant_term

# Round the calculated ferrite level to the nearest 10
rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10

# Print the full equation with the numbers used
print(f"To estimate the ferrite level, we use the formula:")
print(f"Ferrite Level = ({cr_factor} * Cr_eq) - ({ni_factor} * Ni_eq) - {constant_term}")
print(f"Substituting the given values:")
print(f"Ferrite Level = ({cr_factor} * {cr_eq}) - ({ni_factor} * {ni_eq}) - {constant_term}")

# Print the step-by-step calculation result
print(f"Result of calculation: {ferrite_level:.2f}")

# Print the final answer rounded to the nearest 10
print(f"The approximate ferrite level rounded to the nearest 10 is: {rounded_ferrite_level}")