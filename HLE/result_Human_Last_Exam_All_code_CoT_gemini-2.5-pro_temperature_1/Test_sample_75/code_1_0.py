import math

# Step 1 & 2: Define the start and end years of the Second Great War.
start_year = 2504
end_year = 2506

# Step 3: Calculate the duration.
duration = end_year - start_year

# Step 4: Round the duration up to the nearest whole number.
rounded_duration = math.ceil(duration)

# Step 5: Display the equation and the final result.
print(f"The Second Great War in StarCraft lore lasted from {start_year} to {end_year}.")
print("Calculation of the duration (in years):")
print(f"{end_year} - {start_year} = {duration}")
print(f"Rounded up, the war lasted for {rounded_duration} years.")