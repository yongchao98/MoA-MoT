import math

# The limiting factor for the date of the equation is the use of the '@' symbol
# as a mathematical variable. This became plausible after its inclusion in the
# first official ASCII standard.
year = 1963

# The problem requires rounding this year to the nearest 10.
# The numbers involved in this final calculation are the year and 10.
divisor = 10
multiplier = 10

# Perform the rounding calculation.
# We divide the year by 10, round it to the nearest whole number, and then multiply by 10.
rounded_year = round(year / divisor) * multiplier

# Print the final equation with each number used in the calculation, as requested.
print("The final calculation is rounding the earliest plausible year to the nearest 10.")
print(f"round({year} / {divisor}) * {multiplier} = {rounded_year}")