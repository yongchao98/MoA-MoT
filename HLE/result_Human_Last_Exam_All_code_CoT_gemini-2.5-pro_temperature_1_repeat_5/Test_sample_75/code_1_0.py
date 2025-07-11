import math

# The Second Great War begins with the Brood War and ends after the Legacy of the Void epilogue.
start_year = 2500
end_year = 2506

# Calculate the duration
duration = end_year - start_year

# Round up the duration as requested
final_duration = math.ceil(duration)

# Print the explanation and the final equation
print("The Second Great War in StarCraft lore began in the year 2500 and concluded in 2506.")
print(f"To find the duration, we subtract the start year from the end year.")
print(f"Final Equation: {end_year} - {start_year} = {duration}")
print(f"The total duration rounded up is {final_duration} years.")