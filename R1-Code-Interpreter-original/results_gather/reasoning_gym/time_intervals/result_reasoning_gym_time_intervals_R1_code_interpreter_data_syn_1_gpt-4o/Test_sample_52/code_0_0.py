from datetime import datetime

# Define the start and end times
start_time = datetime(2602, 2, 17, 0, 33)
end_time = datetime(2602, 3, 29, 7, 43)

# Calculate the difference
duration = end_time - start_time

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}")