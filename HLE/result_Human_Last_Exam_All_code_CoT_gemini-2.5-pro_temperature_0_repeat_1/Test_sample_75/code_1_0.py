import math

# According to StarCraft lore, the Second Great War began in 2504
# and concluded in 2508.
start_year = 2504
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# Round up the result (though it's already an integer)
rounded_duration = math.ceil(duration)

# Print the final equation
print(f"The Second Great War started in {start_year} and ended in {end_year}.")
print(f"The duration calculation is: {end_year} - {start_year} = {duration}")
print(f"The duration in years, rounded up, is: {rounded_duration}")