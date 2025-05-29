from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 6, 28)
end_date = datetime(2022, 9, 10)

# Initialize a counter for Sundays
sunday_count = 0

# Iterate over each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Sunday (weekday() == 6)
    if current_date.weekday() == 6:
        sunday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the number of Sundays
print(sunday_count)