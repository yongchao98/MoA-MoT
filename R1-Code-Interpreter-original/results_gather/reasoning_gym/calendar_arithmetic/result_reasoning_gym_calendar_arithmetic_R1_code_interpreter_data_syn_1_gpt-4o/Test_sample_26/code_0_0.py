from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 4, 11)
end_date = datetime(2022, 5, 7)

# Initialize a counter for Saturdays
saturday_count = 0

# Iterate through each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a Saturday (5 represents Saturday)
    if current_date.weekday() == 5:
        saturday_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the total number of Saturdays
print(saturday_count)