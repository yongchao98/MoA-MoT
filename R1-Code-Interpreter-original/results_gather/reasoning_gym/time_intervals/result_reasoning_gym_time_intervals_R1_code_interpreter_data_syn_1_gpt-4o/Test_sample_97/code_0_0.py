from datetime import datetime, timedelta
import pytz

# Define the start and end times with their respective time zones
start_time = datetime(1966, 4, 13, 22, 25, 40, tzinfo=pytz.FixedOffset(-300))
end_time = datetime(1966, 5, 12, 10, 25, 16, tzinfo=pytz.FixedOffset(-540))

# Convert both times to UTC
start_time_utc = start_time.astimezone(pytz.utc)
end_time_utc = end_time.astimezone(pytz.utc)

# Calculate the duration
duration = end_time_utc - start_time_utc

# Extract days, hours, minutes, and seconds
days = duration.days
seconds = duration.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")