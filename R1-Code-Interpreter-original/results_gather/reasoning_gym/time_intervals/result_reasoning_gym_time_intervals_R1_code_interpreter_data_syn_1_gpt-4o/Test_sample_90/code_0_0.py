from datetime import datetime
from dateutil import tz

# Define the start and end times with their respective time zones
start_time_str = "2001-11-23 23:23:48 +0900"
end_time_str = "2001-12-22 16:04:33 +1100"

# Parse the datetime strings into datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S %z")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S %z")

# Convert both times to UTC
start_time_utc = start_time.astimezone(tz.UTC)
end_time_utc = end_time.astimezone(tz.UTC)

# Calculate the difference
duration = end_time_utc - start_time_utc

# Extract days, hours, minutes, and seconds
days = duration.days
seconds = duration.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")