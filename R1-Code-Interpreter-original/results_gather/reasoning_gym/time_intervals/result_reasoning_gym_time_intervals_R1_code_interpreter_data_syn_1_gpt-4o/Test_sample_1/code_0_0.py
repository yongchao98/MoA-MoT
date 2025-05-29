from datetime import datetime, timedelta

# Define the start and end times with their respective time zones
start_time_str = "2298-01-14 15:41 -0700"
end_time_str = "2298-03-24 10:53 +0800"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M %z")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M %z")

# Convert both times to UTC
start_time_utc = start_time.astimezone(tz=None)
end_time_utc = end_time.astimezone(tz=None)

# Calculate the difference
duration = end_time_utc - start_time_utc

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result in the required format
print(f"{days} days, {hours:02}:{minutes:02}")