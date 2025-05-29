from datetime import datetime

# Define the departure and arrival dates
departure_date = datetime(2761, 6, 29)
arrival_date = datetime(2761, 7, 30)

# Calculate the difference in days
flight_duration = (arrival_date - departure_date).days

# Print the result
print(flight_duration)