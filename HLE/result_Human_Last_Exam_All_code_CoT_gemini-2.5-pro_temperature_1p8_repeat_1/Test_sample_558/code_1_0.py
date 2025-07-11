import sys
# This script determines the bass singer from the opera riddle.

# Step 1: Identify the opera based on the clues.
opera_title = "Linda di Chamounix"
performer_1 = "Marietta Alboni"
location_1 = "La Scala"
year_1 = 1843

# Step 2: Identify the second production and the time gap.
performer_2 = "Enrico Caruso"
location_2 = "Met Opera"
caruso_last_performance_year = 1905
time_gap_years = 70

# Calculate the target year for the revival.
revival_after_year = caruso_last_performance_year + time_gap_years

# Step 3: Identify the specific revival performance in New York City.
revival_year = 1981
revival_company = "Opera Orchestra of New York"
revival_location = "Carnegie Hall, New York City"

# Step 4: Identify the bass role and the singer in that production.
bass_role = "Il Prefetto"
bass_singer = "Paul Plishka"

# Print the step-by-step reasoning.
print(f"The opera sung by {performer_1} in the title role at {location_1} in {year_1} and also featuring {performer_2} at the {location_2} is '{opera_title}'.")
print(f"{performer_2}'s era in the production was around {caruso_last_performance_year}.")
print(f"The prompt specifies a revival 'more than {time_gap_years} years' later, which means after the year {caruso_last_performance_year} + {time_gap_years} = {revival_after_year}.")
print(f"This specific revival occurred in {revival_year} by the {revival_company} at {revival_location}.")
print(f"In that production, the bass role of {bass_role} was sung by:")
print(f"{bass_singer}")
