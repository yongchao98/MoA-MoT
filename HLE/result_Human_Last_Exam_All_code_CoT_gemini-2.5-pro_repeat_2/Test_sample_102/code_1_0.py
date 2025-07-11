import math

# Define the given nickel and chromium equivalent values
ni_eq = 29
cr_eq = 39

# Define the constants from the Schoefer empirical formula
c1 = 3.34
c2 = 2.46
c3 = 28.6

# Calculate the ferrite percentage using the formula:
# Ferrite % = (c1 * Cr_eq) - (c2 * Ni_eq) - c3
ferrite_percentage = (c1 * cr_eq) - (c2 * ni_eq) - c3

# Round the calculated percentage to the nearest 10
# The second argument of round(), -1, specifies rounding to the nearest 10.
rounded_ferrite = int(round(ferrite_percentage, -1))

# Print the full equation with the final rounded result
print(f"{c1} * {cr_eq} - {c2} * {ni_eq} - {c3} = {rounded_ferrite}")