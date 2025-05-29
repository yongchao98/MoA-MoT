from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 8, 21)
end_date = datetime(2022, 10, 3)

# Initialize a counter for Fridays
friday_count = 0

# Iterate through each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Friday (weekday() returns 4 for Friday)
    if current_date.weekday() == 4:
        friday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the number of Fridays
print(friday_count)