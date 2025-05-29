from datetime import datetime

# Define the departure and arrival times
departure_time = "2139-08-07 02:38"
arrival_time = "2139-09-17 22:38"

# Convert the strings to datetime objects
departure_dt = datetime.strptime(departure_time, "%Y-%m-%d %H:%M")
arrival_dt = datetime.strptime(arrival_time, "%Y-%m-%d %H:%M")

# Calculate the difference
flight_duration = arrival_dt - departure_dt

# Extract days, hours, and minutes
days = flight_duration.days
hours, remainder = divmod(flight_duration.seconds, 3600)
minutes = remainder // 60

# Print the result in the format D days, HH:MM
print(f"{days} days, {hours:02}:{minutes:02}")