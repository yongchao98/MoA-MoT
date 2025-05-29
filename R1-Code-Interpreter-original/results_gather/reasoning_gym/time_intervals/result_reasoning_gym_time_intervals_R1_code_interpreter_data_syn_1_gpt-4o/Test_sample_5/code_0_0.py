from datetime import datetime

# Define the start and end timestamps
start_time = "2805-10-28 11:28:33"
end_time = "2805-11-11 04:06:34"

# Parse the timestamps into datetime objects
start_dt = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
end_dt = datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")

# Calculate the difference
duration = end_dt - start_dt

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, _ = divmod(remainder, 60)

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")