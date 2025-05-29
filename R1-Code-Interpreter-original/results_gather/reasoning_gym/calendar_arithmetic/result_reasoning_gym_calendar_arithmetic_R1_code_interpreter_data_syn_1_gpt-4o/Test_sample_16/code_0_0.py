from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 6, 11)
end_date = datetime(2022, 6, 27)

# Initialize a counter for Fridays
friday_count = 0

# Iterate over each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Friday (weekday() == 4)
    if current_date.weekday() == 4:
        friday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the count of Fridays
print(friday_count)