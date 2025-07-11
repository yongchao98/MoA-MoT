import math

# According to StarCraft lore, the Second Great War started in the year 2500
# and concluded at the end of the Brood War in 2504.
start_year = 2500
end_year = 2504

# Calculate the duration by subtracting the start year from the end year.
duration = end_year - start_year

# The result is an integer, but we use math.ceil() to follow the "rounded up" instruction.
final_duration = math.ceil(duration)

# Print out the explanation and the final calculation.
print(f"The Second Great War in StarCraft Lore began in the year {start_year} and ended in the year {end_year}.")
print(f"The duration is calculated as: {end_year} - {start_year} = {duration}")
print(f"Therefore, the war lasted for {final_duration} years.")
