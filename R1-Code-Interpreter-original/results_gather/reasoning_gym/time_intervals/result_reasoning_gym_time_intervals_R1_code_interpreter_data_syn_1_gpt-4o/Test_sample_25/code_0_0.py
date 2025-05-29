# Define start and end times
start_hour, start_minute = 11, 37
end_hour, end_minute = 11, 58

# Calculate the difference in minutes
duration_minutes = (end_hour * 60 + end_minute) - (start_hour * 60 + start_minute)

# Convert the duration into hours and minutes
duration_hours = duration_minutes // 60
duration_remaining_minutes = duration_minutes % 60

# Print the result in HH:MM format
print(f"{duration_hours:02}:{duration_remaining_minutes:02}")