from datetime import datetime

# Define the departure and arrival dates
departure_date = datetime(2804, 8, 25)
arrival_date = datetime(2804, 11, 11)

# Calculate the difference in days
flight_duration = (arrival_date - departure_date).days

# Print the result
print(flight_duration)