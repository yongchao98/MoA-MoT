from datetime import datetime, timedelta

# Start and end dates
start_date = datetime(2022, 5, 7)
end_date = datetime(2022, 5, 31)

# Initialize the count of Saturdays
saturday_count = 0

# Iterate over each day in the range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Saturday (weekday() == 5 for Saturday)
    if current_date.weekday() == 5:
        saturday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the count of Saturdays
print(saturday_count)