import math

# According to StarCraft lore, the Second Great War began in the year 2504
# with the events of StarCraft II: Wings of Liberty.
start_year = 2504

# The war concluded in the year 2508 with the epilogue of StarCraft II: Legacy of the Void.
end_year = 2508

# Calculate the duration in years.
duration = end_year - start_year

# The problem asks to round up, but since the result is an integer, rounding doesn't change the value.
# We can still use math.ceil to formally adhere to the rounding instruction.
final_duration = math.ceil(duration)

print(f"The Second Great War in StarCraft lore began in {start_year} and concluded in {end_year}.")
print("The calculation for the duration is:")
print(f"{end_year} - {start_year} = {final_duration}")
print(f"So, the war lasted for {final_duration} years.")