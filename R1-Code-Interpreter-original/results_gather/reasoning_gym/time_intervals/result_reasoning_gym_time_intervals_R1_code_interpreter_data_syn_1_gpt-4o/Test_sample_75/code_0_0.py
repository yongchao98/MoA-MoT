# Convert time to seconds since the start of the day
def time_to_seconds(hours, minutes, seconds):
    return hours * 3600 + minutes * 60 + seconds

# Convert seconds back to HH:MM:SS
def seconds_to_hhmmss(total_seconds):
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return hours, minutes, seconds

# Departure and arrival times
departure_time = (19, 37, 49)
arrival_time = (22, 56, 23)

# Convert times to seconds
departure_seconds = time_to_seconds(*departure_time)
arrival_seconds = time_to_seconds(*arrival_time)

# Calculate duration in seconds
duration_seconds = arrival_seconds - departure_seconds

# Convert duration back to HH:MM:SS
duration_hhmmss = seconds_to_hhmmss(duration_seconds)

# Print the result
print(f"{duration_hhmmss[0]:02}:{duration_hhmmss[1]:02}:{duration_hhmmss[2]:02}")