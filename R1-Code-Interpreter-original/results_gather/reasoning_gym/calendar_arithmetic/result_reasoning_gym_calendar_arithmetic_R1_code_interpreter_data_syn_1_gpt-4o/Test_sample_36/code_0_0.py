from datetime import datetime, timedelta

# Define the start and end dates
start_date = datetime(2022, 8, 27)
end_date = datetime(2022, 11, 23)

# Initialize the count of business days
business_days_count = 0

# Iterate over each day in the date range
current_date = start_date
while current_date <= end_date:
    # Check if the current day is a weekday (Monday=0, ..., Sunday=6)
    if current_date.weekday() < 5:  # Monday to Friday are 0 to 4
        business_days_count += 1
    # Move to the next day
    current_date += timedelta(days=1)

# Output the count of business days
print(business_days_count)