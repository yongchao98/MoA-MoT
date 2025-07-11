import math

# The Second Great War in StarCraft lore starts with the Brood War and ends with Wings of Liberty.
start_year = 2500
end_year = 2504

# Calculate the duration
duration = end_year - start_year

# Round up the result as requested
final_duration = math.ceil(duration)

print("Calculating the duration of the Second Great War in StarCraft lore.")
print(f"The war began in {start_year} and ended in {end_year}.")
print("The calculation is:")
print(f"{end_year} - {start_year} = {duration}")
print(f"The total duration rounded up is {final_duration} years.")