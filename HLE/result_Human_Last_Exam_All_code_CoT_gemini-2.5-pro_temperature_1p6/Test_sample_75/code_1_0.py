import math

# The Second Great War begins in the year 2504 with the start of the Wings of Liberty campaign.
start_year = 2504

# The war concludes in the year 2508 with the end of the Legacy of the Void epilogue.
end_year = 2508

# Calculate the duration by subtracting the start year from the end year.
duration = end_year - start_year

# The request asks to round up the result.
# We will use math.ceil to perform the rounding.
final_duration = math.ceil(duration)

# Print the final equation showing the calculation.
print("The Second Great War started in 2504 and ended in 2508.")
print(f"The total duration is calculated as: {end_year} - {start_year} = {final_duration}")
print(f"So, the war lasted for {final_duration} years.")