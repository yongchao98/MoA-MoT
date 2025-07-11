# The earliest plausible year is determined by the release of the TeX typesetting system,
# which produced the font and layout seen in the image.
year = 1978

# The problem asks to round the year to the nearest 10 years.
# We can do this by dividing by 10, rounding to the nearest integer, and then multiplying by 10.
rounded_year = round(year / 10) * 10

print(f"The earliest plausible year is {year}.")
print(f"Rounding {year} to the nearest 10 years gives {rounded_year}.")