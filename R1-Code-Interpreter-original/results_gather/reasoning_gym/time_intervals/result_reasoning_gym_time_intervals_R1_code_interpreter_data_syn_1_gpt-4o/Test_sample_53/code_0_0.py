# Define the start and end times
start_time = "22:39:52"
end_time = "22:54:30"

# Convert the start time to seconds
start_hours, start_minutes, start_seconds = map(int, start_time.split(':'))
start_total_seconds = start_hours * 3600 + start_minutes * 60 + start_seconds

# Convert the end time to seconds
end_hours, end_minutes, end_seconds = map(int, end_time.split(':'))
end_total_seconds = end_hours * 3600 + end_minutes * 60 + end_seconds

# Calculate the duration in seconds
duration_seconds = end_total_seconds - start_total_seconds

# Convert the duration back to hours, minutes, and seconds
duration_hours = duration_seconds // 3600
duration_seconds %= 3600
duration_minutes = duration_seconds // 60
duration_seconds %= 60

# Print the result in HH:MM:SS format
print(f"{duration_hours:02}:{duration_minutes:02}:{duration_seconds:02}")