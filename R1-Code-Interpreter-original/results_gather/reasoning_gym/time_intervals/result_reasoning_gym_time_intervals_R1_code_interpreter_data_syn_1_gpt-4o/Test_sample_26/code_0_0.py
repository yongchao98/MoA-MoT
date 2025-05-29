from datetime import datetime, timedelta

# Define the start and end times with their respective time zones
start_time_str = "2959-10-25 20:12:02 +0000"
end_time_str = "2959-11-19 12:33:46 -0100"

# Parse the start and end times
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S %z")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S %z")

# Convert end time to UTC
end_time_utc = end_time.astimezone(tz=start_time.tzinfo)

# Calculate the duration
duration = end_time_utc - start_time

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}")