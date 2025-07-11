import math

# The Second Great War in StarCraft lore started in the year 2504.
start_year = 2504

# The war is considered to have ended with the main events of Legacy of the Void in 2506.
end_year = 2506

# The question asks for the duration in years, rounded up.
# The most straightforward interpretation that honors the "rounded up" clause
# is to count the number of calendar years the war was active in.
# For example, a war from late 2001 to early 2003 was active in 2001, 2002, and 2003 (3 years).
# The formula for this is: end_year - start_year + 1
duration = (end_year - start_year) + 1

# We use math.ceil to formally round up, although our calculation method already handles this.
# For example, 2506 - 2504 = 2. But the war took place *during* 2504, 2505, and 2506.
# This is a duration of 3 calendar years.
final_duration = math.ceil(duration)

print(f"The Second Great War began in {start_year} and ended in {end_year}.")
print(f"The calculation for the total number of calendar years the war spanned is: {end_year} - {start_year} + 1 = {duration}")
print(f"Therefore, the duration of the war, rounded up, was {final_duration} years.")