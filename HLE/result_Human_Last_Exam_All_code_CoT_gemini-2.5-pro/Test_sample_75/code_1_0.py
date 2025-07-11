import math

# The Second Great War in StarCraft lore started in 2504 with the events of Wings of Liberty.
start_year = 2504

# The war concluded in 2508 with the epilogue of Legacy of the Void.
end_year = 2508

# Calculate the duration.
duration = end_year - start_year

# The problem asks to round up the result.
# Since the duration is an exact number of years, math.ceil() will not change the value.
final_duration = math.ceil(duration)

print("The duration of the Second Great War is calculated as:")
print(f"{end_year} - {start_year} = {final_duration}")
print(f"The total length of the war was {final_duration} years.")
