import math

# According to StarCraft lore, the Second Great War started in the year 2504
# and concluded in the year 2508.
start_year = 2504
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# Round the duration up to the nearest whole number
rounded_up_duration = math.ceil(duration)

print(f"The Second Great War started in {start_year} and ended in {end_year}.")
print(f"Calculation: {end_year} - {start_year} = {duration}")
print(f"The duration of the war in years, rounded up, is: {rounded_up_duration}")