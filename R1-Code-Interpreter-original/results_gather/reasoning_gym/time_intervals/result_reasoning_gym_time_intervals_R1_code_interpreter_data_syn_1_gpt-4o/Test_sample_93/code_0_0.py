from datetime import datetime

# Define the start and end times
start_time = "07:53:51"
end_time = "18:47:45"

# Convert the times to datetime objects
start_dt = datetime.strptime(start_time, "%H:%M:%S")
end_dt = datetime.strptime(end_time, "%H:%M:%S")

# Calculate the difference in seconds
duration_seconds = (end_dt - start_dt).seconds

# Convert the difference back to hours, minutes, and seconds
hours = duration_seconds // 3600
minutes = (duration_seconds % 3600) // 60
seconds = duration_seconds % 60

# Print the duration in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")