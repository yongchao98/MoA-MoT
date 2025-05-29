from datetime import datetime, timedelta

# Given date
start_date = datetime(2022, 10, 6)

# Subtract 88 days
days_to_subtract = 88
new_date = start_date - timedelta(days=days_to_subtract)

# Get the weekday name
weekday_name = new_date.strftime("%A")

print(weekday_name)