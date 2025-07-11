import math

# According to StarCraft lore, the Second Great War began with the Brood War.
start_year = 2500

# The war concluded with the events of the Legacy of the Void epilogue.
end_year = 2508

# Calculate the duration.
# The user requested the result to be rounded up, so we use math.ceil,
# though it's a whole number in this case.
duration = math.ceil(end_year - start_year)

# Print the final equation and the result.
print(f"The Second Great War started in {start_year} and ended in {end_year}.")
print(f"To find the duration, we calculate: {end_year} - {start_year} = {duration}")
print(f"The total duration of the war was {duration} years.")
