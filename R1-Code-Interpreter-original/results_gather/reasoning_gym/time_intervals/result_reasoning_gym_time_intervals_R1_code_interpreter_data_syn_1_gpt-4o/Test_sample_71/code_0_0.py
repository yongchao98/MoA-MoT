from datetime import datetime, timedelta

# Define the start and end times with their respective time zones
start_time_str = "23 February 2889, 02:58 +0800"
end_time_str = "28 May 2889, 07:24 +0000"

# Convert the start and end times to datetime objects
start_time = datetime.strptime(start_time_str, "%d %B %Y, %H:%M %z")
end_time = datetime.strptime(end_time_str, "%d %B %Y, %H:%M %z")

# Calculate the duration
duration = end_time - start_time

# Extract days, hours, and minutes from the duration
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")