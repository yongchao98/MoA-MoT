# 1. Define input values
cr_eq = 39
ni_eq = 29

# 2. Use an empirical formula derived from the WRC-1992 diagram,
#    which is accurate for the specified high Ni_eq level.
#    The formula is: FN = 2 * Cr_eq - 60
#    where FN is the Ferrite Number, approximately the ferrite percentage.
constant_1 = 2
constant_2 = 60
calculated_fn = constant_1 * cr_eq - constant_2

# 3. Round the calculated ferrite level to the nearest 10.
#    The formula for rounding a number 'x' to the nearest 10 is round(x / 10) * 10.
final_ferrite_level = round(calculated_fn / 10) * 10

# 4. Print the output, showing the numbers in the equation and the final result.
print(f"For a Chromium Equivalent of {cr_eq}% and a Nickel Equivalent of {ni_eq}%, the ferrite level is calculated.")
print(f"Calculation based on the derived formula: {constant_1} * {cr_eq} - {constant_2} = {calculated_fn}")
print(f"The calculated ferrite level is approximately {calculated_fn}%.")
print(f"Rounded to the nearest 10, the final ferrite level is: {final_ferrite_level}")