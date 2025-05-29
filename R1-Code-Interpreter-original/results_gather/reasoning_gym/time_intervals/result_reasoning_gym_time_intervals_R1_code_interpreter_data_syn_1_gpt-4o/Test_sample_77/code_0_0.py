from datetime import datetime

# Define the start and end times
start_time_str = "07:29:28.059"
end_time_str = "17:27:12.162"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%H:%M:%S.%f")
end_time = datetime.strptime(end_time_str, "%H:%M:%S.%f")

# Calculate the difference
time_difference = end_time - start_time

# Extract the total seconds and milliseconds
total_seconds = time_difference.total_seconds()
milliseconds = int((total_seconds - int(total_seconds)) * 1000)

# Calculate hours, minutes, and seconds
hours = int(total_seconds // 3600)
minutes = int((total_seconds % 3600) // 60)
seconds = int(total_seconds % 60)

# Format the result
result = f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}"

print(result)