from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 4, 10)
end_date = datetime(2022, 6, 18)

# Initialize the counter for Wednesdays
wednesday_count = 0

# Iterate over each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Wednesday (weekday() == 2)
    if current_date.weekday() == 2:
        wednesday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the count of Wednesdays
print(wednesday_count)