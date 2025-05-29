from datetime import datetime

# Define the start and end times
start_time = datetime(2983, 1, 22, 23, 6)
end_time = datetime(2983, 3, 10, 15, 58)

# Calculate the difference
time_difference = end_time - start_time

# Extract days, hours, minutes, and seconds
days = time_difference.days
seconds = time_difference.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")