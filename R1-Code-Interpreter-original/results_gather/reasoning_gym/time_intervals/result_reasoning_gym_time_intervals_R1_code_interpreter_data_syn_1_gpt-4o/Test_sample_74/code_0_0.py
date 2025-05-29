from datetime import datetime, timedelta

# Define the start and end times in their respective time zones
start_time_str = "30 March 2697 16:52 -1000"
end_time_str = "07 April 2697 02:18 -0600"

# Convert the strings to datetime objects
start_time = datetime.strptime(start_time_str, "%d %B %Y %H:%M %z")
end_time = datetime.strptime(end_time_str, "%d %B %Y %H:%M %z")

# Calculate the difference
time_difference = end_time - start_time

# Extract days, hours, and minutes
days = time_difference.days
hours, remainder = divmod(time_difference.seconds, 3600)
minutes = remainder // 60

# Print the result in the required format
print(f"{days} days, {hours:02}:{minutes:02}")