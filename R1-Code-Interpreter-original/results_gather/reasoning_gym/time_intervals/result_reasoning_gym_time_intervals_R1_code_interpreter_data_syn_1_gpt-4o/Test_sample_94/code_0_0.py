from datetime import datetime

# Define the start and end datetime strings
start_str = "2255-04-07 20:11"
end_str = "2255-06-07 04:31"

# Convert the strings to datetime objects
start_dt = datetime.strptime(start_str, "%Y-%m-%d %H:%M")
end_dt = datetime.strptime(end_str, "%Y-%m-%d %H:%M")

# Calculate the difference
time_difference = end_dt - start_dt

# Extract days, seconds, and then calculate hours, minutes, and seconds
days = time_difference.days
seconds = time_difference.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")