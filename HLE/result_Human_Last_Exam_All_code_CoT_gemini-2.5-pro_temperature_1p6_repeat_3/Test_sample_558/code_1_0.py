# This script solves the riddle by identifying the opera, the specific production, and the singer.

# Define key information from the clues.
year_caruso_last_perf = 1905
years_later_gap = 70
target_year = year_caruso_last_perf + years_later_gap

# Step 1: Identify the opera based on the clues.
# The clues point to Donizetti's "Lucrezia Borgia".
# - Marietta Alboni sang a major role (Maffio Orsini) in a revival at La Scala around 1843.
# - Enrico Caruso sang in it at the Metropolitan Opera.
opera_name = "Lucrezia Borgia"
bass_role = "Don Alfonso d'Este"

# Step 2: Identify the specific NYC production.
# The prompt specifies a revival staged more than 70 years after Caruso's last performance in 1905.
production_company = "New York City Opera"
production_year = 1975
bass_singer = "Richard T. Gill"

# Step 3: Print the step-by-step reasoning.
print("Here is the breakdown of the solution:")
print("--------------------------------------")
print(f"The opera is '{opera_name}' by Gaetano Donizetti.")
print(f"The revival took place more than {years_later_gap} years after Enrico Caruso's last performance in the work.")
print("The calculation for the target year is:")
print(f"{year_caruso_last_perf} + {years_later_gap} = {target_year}")
print(f"\nThis corresponds to the {production_company}'s staged production in {production_year}.")
print(f"The principal bass role in this opera is {bass_role}.")
print(f"In that {production_year} production, the role was sung by:")
print(f"\n{bass_singer}")
