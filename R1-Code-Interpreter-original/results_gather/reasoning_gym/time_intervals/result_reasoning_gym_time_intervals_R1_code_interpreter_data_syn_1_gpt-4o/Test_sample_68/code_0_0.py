# Define the departure and arrival times
departure_time = (14, 29, 49)  # (hours, minutes, seconds)
arrival_time = (17, 40, 42)    # (hours, minutes, seconds)

# Convert times to total seconds from the start of the day
departure_seconds = departure_time[0] * 3600 + departure_time[1] * 60 + departure_time[2]
arrival_seconds = arrival_time[0] * 3600 + arrival_time[1] * 60 + arrival_time[2]

# Calculate the difference in seconds
duration_seconds = arrival_seconds - departure_seconds

# Convert the duration back to hours, minutes, and seconds
hours = duration_seconds // 3600
minutes = (duration_seconds % 3600) // 60
seconds = duration_seconds % 60

# Print the result in HH:MM:SS format
print(f"{hours:02}:{minutes:02}:{seconds:02}")