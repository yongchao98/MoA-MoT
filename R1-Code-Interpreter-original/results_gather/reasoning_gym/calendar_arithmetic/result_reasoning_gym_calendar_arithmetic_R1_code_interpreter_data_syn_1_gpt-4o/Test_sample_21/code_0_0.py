from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 4, 21)
end_date = datetime(2022, 6, 18)

# Find the first Friday on or after the start date
current_date = start_date
while current_date.weekday() != 4:  # 4 corresponds to Friday
    current_date += timedelta(days=1)

# Count the number of Fridays until the end date
friday_count = 0
while current_date <= end_date:
    friday_count += 1
    current_date += timedelta(days=7)

print(friday_count)