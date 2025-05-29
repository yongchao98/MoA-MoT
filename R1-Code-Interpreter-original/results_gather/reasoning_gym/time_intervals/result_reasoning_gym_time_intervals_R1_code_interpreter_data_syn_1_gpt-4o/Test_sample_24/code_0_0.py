from datetime import datetime

# Define the start and end times
start_time = datetime.strptime("2837-06-30 16:19", "%Y-%m-%d %H:%M")
end_time = datetime.strptime("2837-09-08 07:38", "%Y-%m-%d %H:%M")

# Calculate the difference
time_difference = end_time - start_time

# Extract days, hours, minutes, and seconds
days = time_difference.days
seconds = time_difference.seconds
hours = seconds // 3600
minutes = (seconds % 3600) // 60
seconds = seconds % 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}:{seconds:02}")