import math

# A dictionary mapping symbols to their year of introduction or widespread availability for typesetting.
symbol_years = {
    "Equals sign '='": 1557,
    "Plus-minus sign '±'": 1631,
    "Superscript exponent 's²'": 1637,
    "Newton's dot notation 'ṡ'": 1665,
    "Leibniz's derivative notation 'd/ds'": 1684,
    "At sign '@' on a commercial typewriter": 1896
}

# Find the latest year, which is the limiting factor
first_possible_year = 0
for symbol, year in symbol_years.items():
    if year > first_possible_year:
        first_possible_year = year

print(f"The first possible year, based on the introduction of the '@' symbol to typewriters, is {first_possible_year}.")

# Round the year to the nearest 10
# We can do this by dividing by 10, rounding to the nearest integer, and then multiplying by 10.
decade_divisor = 10
rounded_year = round(first_possible_year / decade_divisor) * decade_divisor

print(f"To round this to the nearest 10 years, we perform the calculation:")
print(f"rounded_year = round({first_possible_year} / {decade_divisor}) * {decade_divisor}")
print(f"rounded_year = round({first_possible_year / decade_divisor}) * {decade_divisor}")
print(f"rounded_year = {round(first_possible_year / decade_divisor)} * {decade_divisor}")
print(f"The resulting year is: {rounded_year}")