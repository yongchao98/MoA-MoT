from datetime import datetime

# Define the start and end times
start_time = datetime.strptime("09:44", "%H:%M")
end_time = datetime.strptime("12:22", "%H:%M")

# Calculate the difference
duration = end_time - start_time

# Format the duration in HH:MM
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result
print(f"{hours:02}:{minutes:02}")