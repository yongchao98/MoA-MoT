from datetime import datetime

# Define the start and end times
start_time_str = "02:27:35.584"
end_time_str = "18:58:39.935"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%H:%M:%S.%f")
end_time = datetime.strptime(end_time_str, "%H:%M:%S.%f")

# Calculate the duration
duration = end_time - start_time

# Print the duration in the desired format
print(str(duration))