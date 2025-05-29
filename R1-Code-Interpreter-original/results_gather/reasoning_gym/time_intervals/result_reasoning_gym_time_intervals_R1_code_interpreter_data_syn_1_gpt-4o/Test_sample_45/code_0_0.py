from datetime import datetime

# Define the time format
time_format = "%H:%M:%S.%f"

# Parse the timestamps
time1 = datetime.strptime("10:12:15.842", time_format)
time2 = datetime.strptime("13:43:13.804", time_format)

# Calculate the difference
time_difference = time2 - time1

# Extract the total seconds and convert to milliseconds
total_seconds = time_difference.total_seconds()
milliseconds = int((total_seconds - int(total_seconds)) * 1000)

# Format the result as HH:MM:SS.mmm
hours, remainder = divmod(int(total_seconds), 3600)
minutes, seconds = divmod(remainder, 60)

# Print the result
print(f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}")