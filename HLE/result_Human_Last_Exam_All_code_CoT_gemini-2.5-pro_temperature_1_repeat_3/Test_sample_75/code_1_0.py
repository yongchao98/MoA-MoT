import math

# The Second Great War in StarCraft lore is the primary conflict in the StarCraft II trilogy.
# According to the official timeline, the war began in the year 2504 with the events of Wings of Liberty.
start_year = 2504

# The war concluded with the defeat of the dark god Amon in Legacy of the Void, in the year 2506.
end_year = 2506

# A simple subtraction of the years (2506 - 2504) gives a result of 2.
# However, the conflict spanned across three separate calendar years: 2504, 2505, and 2506.
# Since the problem requires rounding up the duration, we calculate the total span of years the war was active in.
# This is calculated by finding the difference and adding 1 to inclusively count the start year.
duration = (end_year - start_year) + 1

print("The Second Great War in StarCraft lore started in {} and ended in {}.".format(start_year, end_year))
print("To find the war's duration in years (rounded up), we calculate the span of years the conflict was active:")
print("{} - {} + 1 = {}".format(end_year, start_year, duration))
print("\nThe war lasted for {} years.".format(duration))