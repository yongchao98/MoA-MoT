from datetime import datetime, timedelta
import pytz

# Define the start and end times with their respective time zones
start_time_str = "2358-04-20 03:51:19 -0100"
end_time_str = "2358-05-13 11:44:31 -0900"

# Parse the datetime strings into datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S %z")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S %z")

# Convert both times to UTC
start_time_utc = start_time.astimezone(pytz.utc)
end_time_utc = end_time.astimezone(pytz.utc)

# Calculate the difference
duration = end_time_utc - start_time_utc

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}")