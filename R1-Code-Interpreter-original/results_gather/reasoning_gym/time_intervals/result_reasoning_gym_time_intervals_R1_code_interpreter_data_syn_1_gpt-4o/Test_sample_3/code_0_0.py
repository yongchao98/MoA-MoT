from datetime import datetime

# Define the start and end times
start_time_str = "Sun Sep 20 2972 10:55:09"
end_time_str = "Tue Oct 27 2972 02:38:17"

# Define the format of the input strings
time_format = "%a %b %d %Y %H:%M:%S"

# Parse the strings into datetime objects
start_time = datetime.strptime(start_time_str, time_format)
end_time = datetime.strptime(end_time_str, time_format)

# Calculate the difference
duration = end_time - start_time

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, _ = divmod(remainder, 60)

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")