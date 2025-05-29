from datetime import datetime

# Define the departure and arrival times
departure_time = "Mon Mar 02 2229 00:56:24"
arrival_time = "Fri Apr 17 2229 03:17:28"

# Parse the date strings into datetime objects
departure_dt = datetime.strptime(departure_time, "%a %b %d %Y %H:%M:%S")
arrival_dt = datetime.strptime(arrival_time, "%a %b %d %Y %H:%M:%S")

# Calculate the difference between the two datetime objects
flight_duration = arrival_dt - departure_dt

# Extract days, hours, and minutes from the duration
days = flight_duration.days
hours, remainder = divmod(flight_duration.seconds, 3600)
minutes = remainder // 60

# Print the result in the desired format
print(f"{days} days, {hours:02}:{minutes:02}")