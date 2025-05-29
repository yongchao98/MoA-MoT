from datetime import datetime

# Define the two timestamps
time1 = "11:49:32.402"
time2 = "13:25:13.695"

# Define the format
time_format = "%H:%M:%S.%f"

# Parse the timestamps into datetime objects
t1 = datetime.strptime(time1, time_format)
t2 = datetime.strptime(time2, time_format)

# Calculate the difference
duration = t2 - t1

# Get the total seconds and milliseconds
total_seconds = duration.total_seconds()
milliseconds = int((total_seconds - int(total_seconds)) * 1000)

# Convert total seconds to hours, minutes, and seconds
hours = int(total_seconds // 3600)
minutes = int((total_seconds % 3600) // 60)
seconds = int(total_seconds % 60)

# Format the result
result = f"{hours:02}:{minutes:02}:{seconds:02}.{milliseconds:03}"

print(result)