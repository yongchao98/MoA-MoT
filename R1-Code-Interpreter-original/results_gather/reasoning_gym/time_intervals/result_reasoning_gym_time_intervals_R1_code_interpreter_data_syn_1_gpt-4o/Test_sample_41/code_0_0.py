from datetime import datetime

# Define the start and end times
start_time = "03:32:30.472"
end_time = "10:12:35.776"

# Convert the time strings to datetime objects
start_dt = datetime.strptime(start_time, "%H:%M:%S.%f")
end_dt = datetime.strptime(end_time, "%H:%M:%S.%f")

# Calculate the difference
duration = end_dt - start_dt

# Format the duration to include milliseconds
duration_str = str(duration)[:-3]  # Remove the last 3 characters to keep milliseconds

print(duration_str)