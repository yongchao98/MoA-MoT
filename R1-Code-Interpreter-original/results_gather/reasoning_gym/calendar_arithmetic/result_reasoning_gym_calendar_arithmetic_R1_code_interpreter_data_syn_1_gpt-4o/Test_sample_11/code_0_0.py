from datetime import datetime, timedelta

# Given date
start_date = datetime(2022, 5, 23)

# Subtract 98 days
days_to_subtract = timedelta(days=98)
new_date = start_date - days_to_subtract

# Get the full weekday name
weekday_name = new_date.strftime("%A")

print(weekday_name)