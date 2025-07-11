def calculate_pack_years(packs_per_day, years_smoked):
  """
  Calculates the smoking history in pack-years and prints the equation.
  A pack-year is a unit for measuring the amount a person has smoked over a long period of time.
  """
  pack_years = packs_per_day * years_smoked
  print("The patient has a significant smoking history. We can quantify this in pack-years.")
  print(f"The calculation is: {packs_per_day} packs/day * {years_smoked} years = {pack_years} pack-years.")

# Numbers from the clinical case
packs = 2
years = 20

calculate_pack_years(packs, years)