from datetime import datetime

# Define start and end times
start_time = datetime.strptime("12:02", "%H:%M")
end_time = datetime.strptime("22:43", "%H:%M")

# Calculate the difference
duration = end_time - start_time

# Extract hours and minutes from the duration
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")