import math

# According to StarCraft lore, the Second Great War started in 2504 and ended in 2508.
start_year = 2504
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# Round the duration up to the nearest whole number
years_rounded_up = math.ceil(duration)

# Print the equation with each number
print(f"Calculation: {end_year} - {start_year} = {duration}")

# Print the final result
print(f"The Second Great War in StarCraft Lore lasted for {years_rounded_up} years (rounded up).")