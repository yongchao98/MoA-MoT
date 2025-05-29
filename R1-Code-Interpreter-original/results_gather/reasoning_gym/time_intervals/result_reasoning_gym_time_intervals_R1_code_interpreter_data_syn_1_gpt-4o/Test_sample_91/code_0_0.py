# Define start and end times
start_time = "08:05"
end_time = "15:30"

# Convert times to minutes since midnight
start_hours, start_minutes = map(int, start_time.split(':'))
end_hours, end_minutes = map(int, end_time.split(':'))

start_total_minutes = start_hours * 60 + start_minutes
end_total_minutes = end_hours * 60 + end_minutes

# Calculate duration in minutes
duration_minutes = end_total_minutes - start_total_minutes

# Convert duration back to hours and minutes
duration_hours = duration_minutes // 60
duration_remaining_minutes = duration_minutes % 60

# Format the result as HH:MM
duration_formatted = f"{duration_hours:02}:{duration_remaining_minutes:02}"

print(duration_formatted)