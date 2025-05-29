# Convert times to minutes
start_time_minutes = 19 * 60 + 37
end_time_minutes = 23 * 60 + 22

# Calculate the difference in minutes
difference_minutes = end_time_minutes - start_time_minutes

# Convert the difference back to hours and minutes
hours = difference_minutes // 60
minutes = difference_minutes % 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")