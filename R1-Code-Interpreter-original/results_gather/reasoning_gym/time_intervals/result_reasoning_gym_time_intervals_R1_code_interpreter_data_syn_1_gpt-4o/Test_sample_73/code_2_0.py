from datetime import datetime

# Define the start and end times
start_time = "Fri Nov 07 1902 09:33:49"
end_time = "Mon Jan 19 1903 21:44:59"

# Parse the times into datetime objects
start_dt = datetime.strptime(start_time, "%a %b %d %Y %H:%M:%S")
end_dt = datetime.strptime(end_time, "%a %b %d %Y %H:%M:%S")

# Calculate the difference
duration = end_dt - start_dt

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes, _ = divmod(remainder, 60)

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")