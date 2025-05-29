from datetime import datetime, timedelta, timezone

# Define the start and end times with their respective time zones
start_time_str = "10 December 2700, 17:24 -0800"
end_time_str = "12 December 2700, 23:12 +0200"

# Convert the start time to UTC
start_time = datetime.strptime(start_time_str, "%d %B %Y, %H:%M %z")
start_time_utc = start_time.astimezone(timezone.utc)

# Convert the end time to UTC
end_time = datetime.strptime(end_time_str, "%d %B %Y, %H:%M %z")
end_time_utc = end_time.astimezone(timezone.utc)

# Calculate the duration
duration = end_time_utc - start_time_utc

# Format the duration in D days, HH:MM:SS
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")