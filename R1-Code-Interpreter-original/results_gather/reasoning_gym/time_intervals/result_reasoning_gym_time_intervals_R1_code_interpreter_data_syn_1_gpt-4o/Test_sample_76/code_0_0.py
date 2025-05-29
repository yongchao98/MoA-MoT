# Convert time to seconds since the start of the day
def time_to_seconds(hours, minutes, seconds):
    return hours * 3600 + minutes * 60 + seconds

# Departure and arrival times
departure_time = (22, 52, 5)
arrival_time = (23, 52, 12)

# Convert both times to seconds
departure_seconds = time_to_seconds(*departure_time)
arrival_seconds = time_to_seconds(*arrival_time)

# Calculate the duration in seconds
duration_seconds = arrival_seconds - departure_seconds

# Convert the duration back to hours, minutes, and seconds
hours = duration_seconds // 3600
minutes = (duration_seconds % 3600) // 60
seconds = duration_seconds % 60

# Print the result in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")