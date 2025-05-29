# Convert time to seconds since the start of the day
def time_to_seconds(hours, minutes, seconds):
    return hours * 3600 + minutes * 60 + seconds

# Convert seconds back to HH:MM:SS format
def seconds_to_hhmmss(total_seconds):
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"{hours:02}:{minutes:02}:{seconds:02}"

# Departure and arrival times
departure_time = (3, 19, 59)  # 03:19:59
arrival_time = (12, 42, 10)   # 12:42:10

# Convert both times to seconds
departure_seconds = time_to_seconds(*departure_time)
arrival_seconds = time_to_seconds(*arrival_time)

# Calculate the duration in seconds
duration_seconds = arrival_seconds - departure_seconds

# Convert the duration back to HH:MM:SS
flight_duration = seconds_to_hhmmss(duration_seconds)

print(flight_duration)