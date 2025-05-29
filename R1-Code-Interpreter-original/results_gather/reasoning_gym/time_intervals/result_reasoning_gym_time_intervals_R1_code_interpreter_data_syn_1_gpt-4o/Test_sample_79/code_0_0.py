# Define the departure and arrival times
departure_time = (10, 31, 12)  # (hours, minutes, seconds)
arrival_time = (17, 45, 29)    # (hours, minutes, seconds)

# Convert times to seconds since midnight
departure_seconds = departure_time[0] * 3600 + departure_time[1] * 60 + departure_time[2]
arrival_seconds = arrival_time[0] * 3600 + arrival_time[1] * 60 + arrival_time[2]

# Calculate the difference in seconds
duration_seconds = arrival_seconds - departure_seconds

# Convert the difference back to hours, minutes, and seconds
hours = duration_seconds // 3600
minutes = (duration_seconds % 3600) // 60
seconds = duration_seconds % 60

# Print the result in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")