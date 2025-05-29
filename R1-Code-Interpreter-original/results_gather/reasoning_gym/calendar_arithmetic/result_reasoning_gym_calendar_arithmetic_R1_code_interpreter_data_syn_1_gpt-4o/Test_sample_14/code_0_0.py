from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 7, 2)

# Subtract 76 days
new_date = start_date - timedelta(days=76)

# Get the day of the week
day_of_week = new_date.strftime('%A')

print(day_of_week)