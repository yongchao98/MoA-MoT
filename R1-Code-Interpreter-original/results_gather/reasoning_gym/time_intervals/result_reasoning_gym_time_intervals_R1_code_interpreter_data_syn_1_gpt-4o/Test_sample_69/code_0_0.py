from datetime import datetime
from dateutil import tz

# Define the start and end times with their respective time zones
start_time_str = "2570-07-22 01:42:32 -0500"
end_time_str = "2570-08-18 23:37:43 -0300"

# Define the format of the datetime strings
datetime_format = "%Y-%m-%d %H:%M:%S %z"

# Parse the datetime strings into datetime objects
start_time = datetime.strptime(start_time_str, datetime_format)
end_time = datetime.strptime(end_time_str, datetime_format)

# Convert both times to UTC
start_time_utc = start_time.astimezone(tz.UTC)
end_time_utc = end_time.astimezone(tz.UTC)

# Calculate the difference
duration = end_time_utc - start_time_utc

# Extract days, seconds, and format the duration
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")