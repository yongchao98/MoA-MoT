# Convert time to minutes since midnight
def time_to_minutes(hours, minutes):
    return hours * 60 + minutes

# Convert minutes back to hours and minutes
def minutes_to_time(minutes):
    hours = minutes // 60
    minutes = minutes % 60
    return hours, minutes

# Departure and arrival times
departure_time = (18, 11)  # 18:11
arrival_time = (19, 39)    # 19:39

# Convert times to minutes
departure_minutes = time_to_minutes(*departure_time)
arrival_minutes = time_to_minutes(*arrival_time)

# Calculate duration in minutes
duration_minutes = arrival_minutes - departure_minutes

# Convert duration back to hours and minutes
duration_hours, duration_minutes = minutes_to_time(duration_minutes)

# Print the result in HH:MM format
print(f"{duration_hours:02}:{duration_minutes:02}")