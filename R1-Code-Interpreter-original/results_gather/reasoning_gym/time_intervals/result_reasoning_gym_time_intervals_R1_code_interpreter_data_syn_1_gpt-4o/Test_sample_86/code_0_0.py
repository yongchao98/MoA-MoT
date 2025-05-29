from datetime import datetime, timedelta
import pytz

# Define the timestamps and their respective time zones
start_str = "2042-02-17 22:23 -0300"
end_str = "2042-05-26 23:02 -0200"

# Parse the timestamps into datetime objects with time zone information
start = datetime.strptime(start_str, "%Y-%m-%d %H:%M %z")
end = datetime.strptime(end_str, "%Y-%m-%d %H:%M %z")

# Convert both timestamps to UTC
start_utc = start.astimezone(pytz.utc)
end_utc = end.astimezone(pytz.utc)

# Calculate the difference
duration = end_utc - start_utc

# Format the duration in D days, HH:MM:SS
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")