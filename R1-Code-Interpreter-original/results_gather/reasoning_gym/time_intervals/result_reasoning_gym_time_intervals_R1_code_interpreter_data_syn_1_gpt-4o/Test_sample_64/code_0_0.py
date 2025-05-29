from datetime import datetime

# Define the time format
time_format = "%H:%M:%S.%f"

# Parse the timestamps
time1 = datetime.strptime("10:58:43.946", time_format)
time2 = datetime.strptime("18:11:19.156", time_format)

# Calculate the time difference
time_difference = time2 - time1

# Extract hours, minutes, seconds, and milliseconds
hours, remainder = divmod(time_difference.seconds, 3600)
minutes, seconds = divmod(remainder, 60)
milliseconds = time_difference.microseconds // 1000

# Format the result
result = f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}"

print(result)