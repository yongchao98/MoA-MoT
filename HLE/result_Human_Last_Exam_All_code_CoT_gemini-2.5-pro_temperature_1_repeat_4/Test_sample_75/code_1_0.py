import math

# The Second Great War in StarCraft lore started in 2504.
start_year = 2504

# The war is considered to have ended with the defeat of Amon's forces on Aiur in 2506.
end_year = 2506

# Calculate the duration
duration = end_year - start_year

# Round the duration up to the nearest whole number as requested.
# In this case, the result is already an integer.
rounded_duration = math.ceil(duration)

print(f"The Second Great War began in {start_year} and ended in {end_year}.")
print(f"Equation: {end_year} - {start_year} = {duration}")
print(f"The total duration, rounded up, is {rounded_duration} years.")