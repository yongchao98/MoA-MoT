# Convert time to minutes since midnight
def time_to_minutes(hours, minutes):
    return hours * 60 + minutes

# Convert minutes back to HH:MM format
def minutes_to_time(minutes):
    hours = minutes // 60
    minutes = minutes % 60
    return f"{hours:02}:{minutes:02}"

# Departure and arrival times
departure_time = (6, 46)  # 06:46
arrival_time = (10, 55)   # 10:55

# Convert times to minutes
departure_minutes = time_to_minutes(*departure_time)
arrival_minutes = time_to_minutes(*arrival_time)

# Calculate flight duration in minutes
flight_duration_minutes = arrival_minutes - departure_minutes

# Convert duration back to HH:MM format
flight_duration = minutes_to_time(flight_duration_minutes)

print(flight_duration)