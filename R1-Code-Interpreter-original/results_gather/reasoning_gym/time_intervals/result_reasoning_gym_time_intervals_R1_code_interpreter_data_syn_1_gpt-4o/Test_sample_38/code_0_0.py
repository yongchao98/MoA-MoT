# Convert times to minutes since midnight
start_time_minutes = 9 * 60 + 44
end_time_minutes = 17 * 60 + 56

# Calculate the difference in minutes
difference_minutes = end_time_minutes - start_time_minutes

# Convert the difference back to hours and minutes
hours = difference_minutes // 60
minutes = difference_minutes % 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")