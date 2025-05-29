from datetime import datetime

# Define the start and end times
start_time = datetime.strptime('13:51:23', '%H:%M:%S')
end_time = datetime.strptime('14:47:02', '%H:%M:%S')

# Calculate the difference
duration = end_time - start_time

# Print the duration in HH:MM:SS format
print(str(duration))