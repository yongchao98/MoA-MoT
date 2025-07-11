import math

# According to StarCraft lore, the Second Great War began in 2504
# and concluded with the epilogue campaign of Legacy of the Void in 2508.
start_year = 2504
end_year = 2508

# Calculate the duration of the war
duration = end_year - start_year

# The result is a whole number, but we use math.ceil to adhere to the "round up" instruction.
final_duration = math.ceil(duration)

# Output the equation and the final answer
print(f"The Second Great War started in {start_year} and ended in {end_year}.")
print(f"The total duration is calculated as: {end_year} - {start_year} = {duration} years.")
print(f"Rounding up, the war lasted for {final_duration} years.")