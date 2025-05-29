from datetime import datetime, timedelta
import pytz

# Define the start and end times with their respective time zones
start_time_str = "2811-11-24 23:37:42"
end_time_str = "2812-01-23 02:57:01"

# Define the time zones
start_tz = pytz.timezone('Etc/GMT+2')  # UTC-2
end_tz = pytz.timezone('Etc/GMT+6')    # UTC-6

# Parse the datetime strings and localize them to their respective time zones
start_time = start_tz.localize(datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S"))
end_time = end_tz.localize(datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S"))

# Convert both times to UTC
start_time_utc = start_time.astimezone(pytz.utc)
end_time_utc = end_time.astimezone(pytz.utc)

# Calculate the duration
duration = end_time_utc - start_time_utc

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, _ = divmod(remainder, 60)

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")