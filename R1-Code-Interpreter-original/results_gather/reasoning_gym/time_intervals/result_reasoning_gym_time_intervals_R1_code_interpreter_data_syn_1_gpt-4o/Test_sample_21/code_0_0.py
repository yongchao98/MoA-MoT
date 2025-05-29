from datetime import datetime, timedelta

# Define the start and end times with their respective time zones
start_time_str = "2620-09-08 22:30:59 -0100"
end_time_str = "2620-12-09 17:28:10 +1100"

# Convert the time strings to datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S %z")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S %z")

# Convert both times to UTC
start_time_utc = start_time.astimezone(tz=None)
end_time_utc = end_time.astimezone(tz=None)

# Calculate the duration
duration = end_time_utc - start_time_utc

# Extract days, hours, minutes, and seconds from the duration
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the result in the format D days, HH:MM:SS
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")