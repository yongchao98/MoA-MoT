import math

# The Second Great War in StarCraft lore begins in 2504.
start_year = 2504

# The war's final epilogue, "Into the Void," concludes in 2508.
end_year = 2508

# Calculate the duration of the war.
duration = end_year - start_year

# The result must be rounded up as per the request.
# math.ceil rounds a number up to the nearest integer.
final_duration = math.ceil(duration)

print("The Second Great War began in {} and ended in {}.".format(start_year, end_year))
print("Calculating the duration (in years, rounded up):")
# We use int() to display the final number without the .0 from ceiling.
print("{} - {} = {}".format(end_year, start_year, int(final_duration)))
<<<4>>>