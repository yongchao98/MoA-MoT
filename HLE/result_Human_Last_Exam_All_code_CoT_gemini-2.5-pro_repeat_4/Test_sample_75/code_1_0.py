import math

# According to StarCraft lore, the Second Great War began in 2504 and ended in 2506.
start_year = 2504
end_year = 2506

# Calculate the duration in years.
duration = end_year - start_year

# Round the result up to the nearest whole number.
final_duration = math.ceil(duration)

# Print the final equation showing all the numbers.
print(f"The Second Great War lasted from {start_year} to {end_year}.")
print(f"The calculation for the duration is: {end_year} - {start_year} = {final_duration} years.")