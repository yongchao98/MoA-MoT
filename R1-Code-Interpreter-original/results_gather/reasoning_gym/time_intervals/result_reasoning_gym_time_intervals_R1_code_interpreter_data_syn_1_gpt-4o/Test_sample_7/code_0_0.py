from datetime import datetime

# Define the start and end times
start_time = datetime.strptime('08:24:01', '%H:%M:%S')
end_time = datetime.strptime('11:57:05', '%H:%M:%S')

# Calculate the difference
duration = end_time - start_time

# Print the duration in HH:MM:SS format
print(str(duration))