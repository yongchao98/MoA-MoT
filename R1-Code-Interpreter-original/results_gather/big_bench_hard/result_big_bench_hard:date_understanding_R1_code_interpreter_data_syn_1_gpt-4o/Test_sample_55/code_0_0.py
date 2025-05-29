from datetime import datetime, timedelta

# Start from January 1, 2019
start_date = datetime(2019, 1, 1)

# Find the first Monday of 2019
# January 1, 2019 is a Tuesday, so the first Monday is 6 days later
first_monday = start_date + timedelta(days=6)

# Find the date of the day before the first Monday
day_before_monday = first_monday - timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = day_before_monday.strftime("%m/%d/%Y")
print(formatted_date)