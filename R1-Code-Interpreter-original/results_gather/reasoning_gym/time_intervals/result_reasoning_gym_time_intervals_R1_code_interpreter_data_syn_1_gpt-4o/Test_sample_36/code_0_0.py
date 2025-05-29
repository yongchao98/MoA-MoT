from datetime import datetime, timedelta

# Define the start and end times with their respective time zones
start_time_str = "24 November 2770, 18:15 +1100"
end_time_str = "24 January 2771, 03:55 +0900"

# Convert the strings to datetime objects
start_time = datetime.strptime(start_time_str, "%d %B %Y, %H:%M %z")
end_time = datetime.strptime(end_time_str, "%d %B %Y, %H:%M %z")

# Convert both times to UTC
start_time_utc = start_time.astimezone(tz=None)
end_time_utc = end_time.astimezone(tz=None)

# Calculate the difference
time_difference = end_time_utc - start_time_utc

# Extract days, hours, and minutes
days = time_difference.days
hours, remainder = divmod(time_difference.seconds, 3600)
minutes = remainder // 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")