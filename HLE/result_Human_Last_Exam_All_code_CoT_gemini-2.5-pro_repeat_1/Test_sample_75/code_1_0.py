import math

# According to StarCraft lore, the Second Great War began in 2504
# and concluded in 2506.
start_year = 2504
end_year = 2506

# Calculate the duration
duration = end_year - start_year

# Round up the duration as requested
final_duration = math.ceil(duration)

# Print the equation
print(f"End Year - Start Year = Duration")
print(f"{end_year} - {start_year} = {duration}")

# Print the final answer
print(f"\nThe Second Great War lasted for {final_duration} years.")