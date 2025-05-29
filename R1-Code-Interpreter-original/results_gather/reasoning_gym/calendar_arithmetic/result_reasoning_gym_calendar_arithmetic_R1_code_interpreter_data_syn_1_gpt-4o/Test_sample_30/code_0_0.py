from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 3, 30)
end_date = datetime(2022, 4, 7)

# Initialize the count of Sundays
sunday_count = 0

# Iterate through each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Sunday (weekday() == 6)
    if current_date.weekday() == 6:
        sunday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the count of Sundays
print(sunday_count)