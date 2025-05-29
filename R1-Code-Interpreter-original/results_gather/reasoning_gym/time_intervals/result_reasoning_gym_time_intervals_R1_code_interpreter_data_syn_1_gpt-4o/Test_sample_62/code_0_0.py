from datetime import datetime

# Define the departure and arrival times in UTC
departure_utc = datetime.strptime("2302-01-23 10:27", "%Y-%m-%d %H:%M")
arrival_utc = datetime.strptime("2302-01-26 06:36", "%Y-%m-%d %H:%M")

# Calculate the duration
duration = arrival_utc - departure_utc

# Extract days, hours, and minutes
days = duration.days
hours, remainder = divmod(duration.seconds, 3600)
minutes = remainder // 60

# Print the result in the format D days, HH:MM
print(f"{days} days, {hours:02}:{minutes:02}")