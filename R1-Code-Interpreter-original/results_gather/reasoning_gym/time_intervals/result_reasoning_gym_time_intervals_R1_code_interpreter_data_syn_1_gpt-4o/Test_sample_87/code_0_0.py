from datetime import datetime

# Define the start and end times
start_time_str = "05:47:50.551"
end_time_str = "22:45:21.232"

# Define the format for parsing the time strings
time_format = "%H:%M:%S.%f"

# Parse the time strings into datetime objects
start_time = datetime.strptime(start_time_str, time_format)
end_time = datetime.strptime(end_time_str, time_format)

# Calculate the duration
duration = end_time - start_time

# Format the duration as HH:MM:SS.mmm
duration_str = str(duration)

print(duration_str)