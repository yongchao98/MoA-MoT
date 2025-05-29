from datetime import datetime

# Define the start and end times
start_time_str = "06:49:07.114"
end_time_str = "16:36:03.456"

# Define the format for parsing the time strings
time_format = "%H:%M:%S.%f"

# Parse the time strings into datetime objects
start_time = datetime.strptime(start_time_str, time_format)
end_time = datetime.strptime(end_time_str, time_format)

# Calculate the duration
duration = end_time - start_time

# Format the duration into HH:MM:SS.mmm
duration_str = str(duration)

print(duration_str)