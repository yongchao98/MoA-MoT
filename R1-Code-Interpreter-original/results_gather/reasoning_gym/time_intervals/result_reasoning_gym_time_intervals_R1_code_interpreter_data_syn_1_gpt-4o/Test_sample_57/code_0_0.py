from datetime import datetime

# Define start and end times
start_time = datetime.strptime("19:45", "%H:%M")
end_time = datetime.strptime("20:23", "%H:%M")

# Calculate the difference
time_difference = end_time - start_time

# Extract hours and minutes from the time difference
hours, remainder = divmod(time_difference.seconds, 3600)
minutes = remainder // 60

# Print the result in HH:MM format
print(f"{hours:02}:{minutes:02}")