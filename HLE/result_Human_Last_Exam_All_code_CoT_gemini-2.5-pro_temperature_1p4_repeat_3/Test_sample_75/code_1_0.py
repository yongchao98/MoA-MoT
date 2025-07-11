import math

# The Second Great War in StarCraft lore started in 2504.
start_year = 2504

# The war concluded with the defeat of Amon in the epilogue of Legacy of the Void in 2508.
end_year = 2508

# Calculate the duration in years.
duration = end_year - start_year

# The problem asks to round up the result.
# Since the result is a whole number, rounding it up does not change the value.
final_duration = math.ceil(duration)

# Print the equation as requested.
print(f"The duration of the Second Great War is calculated as: {end_year} - {start_year} = {final_duration} years.")