from datetime import datetime
import pytz

# Define the start and end times with their respective time zones
start_time_str = "2609-01-28 00:08:51 -1000"
end_time_str = "2609-01-29 22:02:04 +0800"

# Convert the time strings to datetime objects with the correct time zones
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S %z")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S %z")

# Convert both times to UTC using pytz
start_time_utc = start_time.astimezone(pytz.utc)
end_time_utc = end_time.astimezone(pytz.utc)

# Calculate the duration
duration = end_time_utc - start_time_utc

# Extract days, hours, and minutes from the duration
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Format the result
result = f"{days} days, {hours:02}:{minutes:02}"
print(result)