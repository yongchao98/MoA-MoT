# The Second Great War in StarCraft lore is a prolonged conflict.
# For the purpose of this calculation, we'll consider its start to be the beginning
# of the Brood War, which set the stage for the events of StarCraft II.

start_year = 2500

# The war's main conflict concludes with the defeat of the dark god Amon in the
# Legacy of the Void epilogue.

end_year = 2506

# Calculate the total duration in years.
duration = end_year - start_year

# The request is to round up, but since the result is a whole number,
# rounding up doesn't change it.
rounded_duration = duration

# Print the equation and the final answer.
print(f"The duration of the Second Great War is calculated from the end year minus the start year.")
print(f"Equation: {end_year} - {start_year} = {rounded_duration}")
print(f"The war lasted for {rounded_duration} years.")
