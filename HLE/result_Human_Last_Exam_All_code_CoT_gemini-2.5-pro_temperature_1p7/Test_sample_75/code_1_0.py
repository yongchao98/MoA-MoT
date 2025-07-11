import math

# The Second Great War in StarCraft lore begins with the events of Wings of Liberty
# and concludes with the epilogue of Legacy of the Void.
start_year = 2504
end_year = 2508

# Calculate the duration in years
duration = end_year - start_year

# The request is to round up the result.
# math.ceil() is used to round a number up to the nearest integer.
final_duration = math.ceil(duration)

print("The Second Great War began in {} and concluded in {}.".format(start_year, end_year))
print("Calculating the duration: {} - {} = {} years.".format(end_year, start_year, duration))
print("Rounding up, the war lasted for a total of {} years.".format(final_duration))