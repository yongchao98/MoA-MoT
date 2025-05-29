from datetime import datetime

# Define start and end times
start_time = datetime.strptime('16:35:10', '%H:%M:%S')
end_time = datetime.strptime('18:49:34', '%H:%M:%S')

# Calculate the difference
duration = end_time - start_time

# Print the duration in HH:MM:SS format
print(str(duration))