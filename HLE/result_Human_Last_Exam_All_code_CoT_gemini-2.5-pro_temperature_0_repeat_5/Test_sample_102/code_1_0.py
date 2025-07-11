import math

# Given equivalent values
cr_eq = 39
ni_eq = 29

# Coefficients from Hull's formula for ferrite prediction
coeff_cr = 3.34
coeff_ni = 2.46
constant = 28.6

# Calculate the ferrite level using Hull's formula
# Ferrite Level = (coeff_cr * Cr_eq) - (coeff_ni * Ni_eq) - constant
term1 = coeff_cr * cr_eq
term2 = coeff_ni * ni_eq
ferrite_level = term1 - term2 - constant

# Round the result to the nearest 10
rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10

# Print the equation and the final result
print(f"Using Hull's formula: (A * Cr_eq) - (B * Ni_eq) - C")
print(f"Equation: ({coeff_cr} * {cr_eq}) - ({coeff_ni} * {ni_eq}) - {constant} = {ferrite_level:.2f}")
print(f"\nThe calculated ferrite level is approximately {ferrite_level:.2f}%.")
print(f"Rounded to the nearest 10, the ferrite level is: {rounded_ferrite_level}")
