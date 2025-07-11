import math

# The Second Great War in StarCraft lore starts with the Brood War and ends with the conclusion of Legacy of the Void.
start_year = 2500
end_year = 2506

# Calculate the duration
duration = end_year - start_year

# Round the duration up to the nearest whole number
rounded_duration = math.ceil(duration)

# Print the explanation and the final equation
print("The Second Great War began in the year 2500 and concluded in 2506.")
print(f"Calculation: {end_year} - {start_year} = {duration}")
print(f"The total duration, rounded up, is {rounded_duration} years.")