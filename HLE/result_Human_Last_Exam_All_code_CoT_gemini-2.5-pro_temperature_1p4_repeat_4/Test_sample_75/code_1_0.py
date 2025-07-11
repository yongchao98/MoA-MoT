import math

# Define the start and end years of the Second Great War
start_year = 2500
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# The request is to round up the result. Since it's a whole number,
# ceiling will return the same integer.
final_duration = math.ceil(duration)

# Print the final equation as requested
print(f"The duration of the Second Great War is calculated as:")
print(f"{end_year} - {start_year} = {final_duration}")
print(f"The war lasted for {final_duration} years.")