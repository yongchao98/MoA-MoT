from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 3, 27)
end_date = datetime(2022, 6, 15)

# Find the first Tuesday on or after the start date
current_date = start_date
while current_date.weekday() != 1:  # 1 corresponds to Tuesday
    current_date += timedelta(days=1)

# Count the number of Tuesdays
tuesday_count = 0
while current_date <= end_date:
    tuesday_count += 1
    current_date += timedelta(weeks=1)

print(tuesday_count)