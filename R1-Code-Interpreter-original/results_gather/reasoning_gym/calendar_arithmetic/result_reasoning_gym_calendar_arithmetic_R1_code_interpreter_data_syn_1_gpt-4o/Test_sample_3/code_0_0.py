from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 6, 21)
end_date = datetime(2022, 9, 25)

# Initialize the count of Tuesdays
tuesday_count = 0

# Iterate through each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Tuesday (weekday() == 1 for Tuesday)
    if current_date.weekday() == 1:
        tuesday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the number of Tuesdays
print(tuesday_count)