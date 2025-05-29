# Convert time to minutes since the start of the day
start_time = 20 * 60 + 55  # 20:55
end_time = 22 * 60 + 10    # 22:10

# Calculate the duration in minutes
duration_minutes = end_time - start_time

# Convert the duration back to hours and minutes
hours = duration_minutes // 60
minutes = duration_minutes % 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")