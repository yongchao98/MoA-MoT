import math

def round_to_nearest_ten(n):
  """Rounds a number to the nearest 10."""
  return int(round(n / 10.0)) * 10

# The limiting factor for the equation's publication date is the invention of the TeX
# typesetting system, which made it practical to use non-standard symbols like '@'
# in mathematical formulas. TeX was created in 1978.
first_possible_year = 1978

# Round the year to the nearest 10 years.
rounded_year = round_to_nearest_ten(first_possible_year)

print(f"The earliest possible year is determined by the introduction of the TeX typesetting system in {first_possible_year}.")
print(f"Rounding {first_possible_year} to the nearest 10 years gives: {rounded_year}")
print(f"Final Answer: {rounded_year}")
<<<1980>>>