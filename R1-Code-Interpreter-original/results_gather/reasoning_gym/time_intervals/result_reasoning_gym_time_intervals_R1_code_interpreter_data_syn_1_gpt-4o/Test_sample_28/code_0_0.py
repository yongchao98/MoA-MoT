from datetime import datetime

# Define the start and end times
start_time_str = "09:36:09.235"
end_time_str = "14:26:18.062"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%H:%M:%S.%f")
end_time = datetime.strptime(end_time_str, "%H:%M:%S.%f")

# Calculate the difference
duration = end_time - start_time

# Print the duration in the desired format
print(duration)