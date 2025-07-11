import math

# 1. Define the chromium and nickel equivalent values.
cr_eq = 39
ni_eq = 29

# 2. Calculate the ferrite percentage using the Espy formula,
#    which is a mathematical approximation of the WRC-1992 diagram.
#    Formula: Ferrite % = 8.16 * CrEq - 3.47 * NiEq - 0.057 * NiEq^2 - 108.6
ferrite_percentage = 8.16 * cr_eq - 3.47 * ni_eq - 0.057 * (ni_eq**2) - 108.6

# 3. Round the calculated ferrite percentage to the nearest 10.
rounded_ferrite = round(ferrite_percentage / 10.0) * 10

# 4. Print the full equation with the original numbers and the final rounded result.
print(f"8.16 * {cr_eq} - 3.47 * {ni_eq} - 0.057 * {ni_eq}^2 - 108.6 = {rounded_ferrite}")
