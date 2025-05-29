from datetime import datetime

# Define the departure and arrival dates
departure_date = "03/10/2231"
arrival_date = "03/31/2231"

# Convert the string dates to datetime objects
departure = datetime.strptime(departure_date, "%m/%d/%Y")
arrival = datetime.strptime(arrival_date, "%m/%d/%Y")

# Calculate the difference in days
duration = (arrival - departure).days

# Print the result
print(f"{duration} days")