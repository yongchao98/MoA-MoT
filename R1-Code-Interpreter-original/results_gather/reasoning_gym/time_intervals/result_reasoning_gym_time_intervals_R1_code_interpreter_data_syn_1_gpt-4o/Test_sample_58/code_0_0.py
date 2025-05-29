from datetime import datetime

# Define the start and end times
start_time = datetime.strptime('22:05:30', '%H:%M:%S')
end_time = datetime.strptime('23:35:33', '%H:%M:%S')

# Calculate the difference
duration = end_time - start_time

# Extract hours, minutes, and seconds from the duration
hours, remainder = divmod(duration.seconds, 3600)
minutes, seconds = divmod(remainder, 60)

# Print the result in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")