from datetime import datetime

# Define the start and end times
start_time_str = "17:00:11.426"
end_time_str = "17:32:34.712"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%H:%M:%S.%f")
end_time = datetime.strptime(end_time_str, "%H:%M:%S.%f")

# Calculate the difference
duration = end_time - start_time

# Extract hours, minutes, seconds, and milliseconds from the duration
hours, remainder = divmod(duration.seconds, 3600)
minutes, seconds = divmod(remainder, 60)
milliseconds = duration.microseconds // 1000

# Format the duration as HH:MM:SS.mmm
formatted_duration = f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}"

print(formatted_duration)