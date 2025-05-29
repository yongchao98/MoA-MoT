# Convert times to minutes since midnight
time1_minutes = 6 * 60 + 9
time2_minutes = 11 * 60 + 40

# Calculate the difference in minutes
difference_minutes = time2_minutes - time1_minutes

# Convert the difference back to hours and minutes
hours = difference_minutes // 60
minutes = difference_minutes % 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")