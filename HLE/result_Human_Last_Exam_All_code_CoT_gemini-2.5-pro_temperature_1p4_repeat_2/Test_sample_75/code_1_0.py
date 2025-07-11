import math

# According to StarCraft lore, the Second Great War began in 2504
# and was formally declared over in 2508.
start_year = 2504
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# The request is to round up, though the result is a whole number.
# We will apply math.ceil() for correctness as per the prompt.
final_duration = math.ceil(duration)

# Print the final equation and result
print(f"The Second Great War in StarCraft lore started in the year {start_year} and officially ended in {end_year}.")
print(f"The calculation for the duration is: {end_year} - {start_year} = {duration}")
print(f"The total duration was {final_duration} years.")