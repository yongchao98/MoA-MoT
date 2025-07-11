import math

# According to StarCraft lore, the Second Great War starts in 2504
# and concludes in the epilogue in 2508.
start_year = 2504
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# The problem asks to round up, which doesn't change an integer result.
final_duration = math.ceil(duration)

print(f"The Second Great War in StarCraft lore started in the year {start_year} and ended in the year {end_year}.")
print(f"To find the duration, we calculate the difference.")
print(f"Calculation: {end_year} - {start_year} = {final_duration}")
print(f"The war lasted for {final_duration} years.")