from datetime import datetime

# Define the start and end times
start_time_str = "Fri Nov 07 1902 09:33:49"
end_time_str = "Mon Jan 19 1903 21:44:59"

# Parse the strings into datetime objects
start_time = datetime.strptime(start_time_str, "%a %b %d %Y %H:%M:%S")
end_time = datetime.strptime(end_time_str, "%a %b %d %Y %H:%M:%S")

# Calculate the difference
duration = end_time - start_time

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")