from datetime import datetime, timedelta

# Given date
start_date = datetime(2022, 4, 5)

# Subtract 39 days
previous_date = start_date - timedelta(days=39)

# Get the weekday name
weekday_name = previous_date.strftime('%A')

print(weekday_name)