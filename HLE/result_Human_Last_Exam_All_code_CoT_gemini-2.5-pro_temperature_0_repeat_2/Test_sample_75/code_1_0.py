import math

# The Second Great War in StarCraft lore starts with the events of Wings of Liberty
# and ends with the epilogue of Legacy of the Void.
start_year = 2504
end_year = 2506

# Calculate the duration in years
duration = end_year - start_year

# The question asks to round up the result.
# Since the result is an integer, ceiling won't change it, but we follow the instruction.
final_duration = math.ceil(duration)

print(f"The Second Great War in StarCraft lore began in the year {start_year} and ended in the year {end_year}.")
print(f"The calculation for the duration is: {end_year} - {start_year} = {duration}")
print(f"Rounding up to the nearest whole year, the war lasted for {final_duration} years.")