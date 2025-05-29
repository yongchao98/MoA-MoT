from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 8, 6)
end_date = datetime(2022, 9, 9)

# Find the first Sunday on or after the start date
current_date = start_date
while current_date.weekday() != 6:  # 6 corresponds to Sunday
    current_date += timedelta(days=1)

# Count the number of Sundays until the end date
sunday_count = 0
while current_date <= end_date:
    sunday_count += 1
    current_date += timedelta(days=7)

print(sunday_count)