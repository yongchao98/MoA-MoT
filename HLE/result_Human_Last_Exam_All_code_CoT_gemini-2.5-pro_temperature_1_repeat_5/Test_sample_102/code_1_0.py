import math

# Define the input equivalents
cr_eq = 39
ni_eq = 29

# Constants from the WRC-1992 approximation formula
factor_cr = 3.34
factor_ni = 2.46
constant = 28.6

# Calculate the ferrite level
ferrite_level = (factor_cr * cr_eq) - (factor_ni * ni_eq) - constant

# Round the result to the nearest 10
# The second argument to round(), -1, rounds to the nearest 10
rounded_ferrite_level = int(round(ferrite_level, -1))

# Print the full equation and the final answer
print(f"The calculation is based on the formula: (3.34 * Chromium Equivalent) - (2.46 * Nickel Equivalent) - 28.6")
print(f"Plugging in the values: ({factor_cr} * {cr_eq}) - ({factor_ni} * {ni_eq}) - {constant} = {ferrite_level:.2f}")
print(f"The result {ferrite_level:.2f} rounded to the nearest 10 is: {rounded_ferrite_level}")