from datetime import datetime, timedelta
import pytz

# Define the timestamps and their respective time zones
timestamp1 = "2241-10-17 14:19 -0700"
timestamp2 = "2242-01-11 03:17 +0400"

# Parse the timestamps into datetime objects with timezone information
dt1 = datetime.strptime(timestamp1, "%Y-%m-%d %H:%M %z")
dt2 = datetime.strptime(timestamp2, "%Y-%m-%d %H:%M %z")

# Convert both timestamps to UTC
dt1_utc = dt1.astimezone(pytz.utc)
dt2_utc = dt2.astimezone(pytz.utc)

# Calculate the difference between the two UTC timestamps
time_difference = dt2_utc - dt1_utc

# Extract days, hours, and minutes from the time difference
days = time_difference.days
hours, remainder = divmod(time_difference.seconds, 3600)
minutes = remainder // 60

# Print the result in the format D days, HH:MM
print(f"{days} days, {hours:02}:{minutes:02}")