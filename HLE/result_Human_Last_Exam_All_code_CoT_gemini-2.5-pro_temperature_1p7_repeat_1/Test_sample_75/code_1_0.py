import math

# The Second Great War in StarCraft lore began with the Brood War in the year 2500.
start_year = 2500

# The war concluded with the end of the conflict against Amon in Legacy of the Void in the year 2506.
end_year = 2506

# Calculate the duration by subtracting the start year from the end year.
# The result is rounded up as requested, although it's a whole number here.
duration = math.ceil(end_year - start_year)

# Print the final result including the equation used for the calculation.
print("The Second Great War in StarCraft lore is generally considered to have started in the year 2500 and ended in 2506.")
print(f"The calculation for the total duration is: {end_year} - {start_year} = {int(duration)}")
print(f"Therefore, the war lasted for {int(duration)} years.")
