import math

# Define the chromium and nickel equivalent values from the user's request.
cr_eq = 39.0
ni_eq = 29.0

# An established empirical formula (Hull's formula) is used to approximate the
# Ferrite Number (FN) from the chromium and nickel equivalents. The FN is
# considered a close approximation of the ferrite percentage for these alloys.
# The formula is: FN = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6
ferrite_level = 3.34 * cr_eq - 2.46 * ni_eq - 28.6

# Round the calculated ferrite level to the nearest 10 as requested.
# We achieve this by dividing by 10, rounding to the nearest whole number,
# and then multiplying by 10.
rounded_ferrite_level = int(round(ferrite_level / 10.0) * 10)

# Print the explanation, the formula with the numbers plugged in, and the final result.
print("To estimate the ferrite level, we use the empirical Hull's formula:")
print("Ferrite Level = 3.34 * (Chromium Equivalent) - 2.46 * (Nickel Equivalent) - 28.6")
print("\nPlugging in the values gives the following equation:")
print(f"Ferrite Level = 3.34 * {cr_eq} - 2.46 * {ni_eq} - 28.6")
print(f"The calculated ferrite level is approximately: {ferrite_level:.2f}")
print("\nRounding to the nearest 10, the final approximate ferrite level is:")
print(rounded_ferrite_level)