# Convert times to minutes since midnight
start_time_minutes = 6 * 60 + 37
end_time_minutes = 13 * 60 + 52

# Calculate the difference in minutes
elapsed_minutes = end_time_minutes - start_time_minutes

# Convert the difference back to hours and minutes
hours = elapsed_minutes // 60
minutes = elapsed_minutes % 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")