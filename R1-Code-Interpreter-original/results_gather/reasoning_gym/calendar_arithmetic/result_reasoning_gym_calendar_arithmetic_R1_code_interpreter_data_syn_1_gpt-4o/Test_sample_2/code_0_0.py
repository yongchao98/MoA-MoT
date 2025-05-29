from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 7, 2)
end_date = datetime(2022, 9, 14)

# Initialize a counter for Mondays
monday_count = 0

# Iterate through each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Monday (0 = Monday, 6 = Sunday)
    if current_date.weekday() == 0:
        monday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the number of Mondays
print(monday_count)