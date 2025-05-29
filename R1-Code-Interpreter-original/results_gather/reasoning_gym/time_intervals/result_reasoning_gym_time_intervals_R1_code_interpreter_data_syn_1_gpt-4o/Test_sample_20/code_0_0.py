from datetime import datetime

# Define the start and end times
start_time_str = "01:56:46.495"
end_time_str = "20:36:33.370"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%H:%M:%S.%f")
end_time = datetime.strptime(end_time_str, "%H:%M:%S.%f")

# Calculate the difference in time
duration = end_time - start_time

# Format the duration as HH:MM:SS.mmm
duration_str = str(duration)

print(duration_str)