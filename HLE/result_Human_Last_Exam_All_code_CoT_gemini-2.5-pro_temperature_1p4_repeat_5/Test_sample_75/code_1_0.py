import math

# The Second Great War in StarCraft lore starts in 2504 and ends in 2508.
start_year = 2504
end_year = 2508

# Calculate the duration
duration = end_year - start_year

# The problem asks to round up, although the result is already a whole number.
# We will use math.ceil to formally meet this requirement.
rounded_duration = math.ceil(duration)

print("The Second Great War in StarCraft lore began in the year {} and ended in the year {}.".format(start_year, end_year))
print("Calculating the duration:")
print("{} - {} = {}".format(end_year, start_year, duration))
print("The total duration is {} years.".format(rounded_duration))