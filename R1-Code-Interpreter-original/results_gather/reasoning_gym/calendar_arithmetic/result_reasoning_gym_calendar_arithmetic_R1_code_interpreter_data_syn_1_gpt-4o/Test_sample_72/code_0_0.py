from datetime import datetime, timedelta

# Given date
start_date = datetime(2022, 4, 24)

# Subtract 99 days
days_to_subtract = timedelta(days=99)
new_date = start_date - days_to_subtract

# Get the full weekday name
weekday_name = new_date.strftime("%A")

print(weekday_name)