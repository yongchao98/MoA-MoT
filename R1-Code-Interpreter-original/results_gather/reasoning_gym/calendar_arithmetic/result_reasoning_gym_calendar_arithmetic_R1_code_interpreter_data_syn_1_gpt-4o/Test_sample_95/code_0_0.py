from datetime import datetime, timedelta

# Given date
given_date = datetime(2022, 10, 19)

# Subtract 73 days
days_ago = given_date - timedelta(days=73)

# Get the weekday name
weekday_name = days_ago.strftime('%A')

print(weekday_name)