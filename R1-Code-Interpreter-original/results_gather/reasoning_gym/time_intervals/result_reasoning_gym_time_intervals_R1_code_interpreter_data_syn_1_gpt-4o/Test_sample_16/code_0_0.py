from datetime import datetime

# Define the departure and arrival times
departure_time = "Tue May 30 2699 22:34:19"
arrival_time = "Thu Aug 03 2699 05:35:19"

# Convert the strings to datetime objects
departure_dt = datetime.strptime(departure_time, "%a %b %d %Y %H:%M:%S")
arrival_dt = datetime.strptime(arrival_time, "%a %b %d %Y %H:%M:%S")

# Calculate the difference
flight_duration = arrival_dt - departure_dt

# Extract days, hours, and minutes
days = flight_duration.days
hours, remainder = divmod(flight_duration.seconds, 3600)
minutes = remainder // 60

# Print the result
print(f"{days} days, {hours:02}:{minutes:02}")