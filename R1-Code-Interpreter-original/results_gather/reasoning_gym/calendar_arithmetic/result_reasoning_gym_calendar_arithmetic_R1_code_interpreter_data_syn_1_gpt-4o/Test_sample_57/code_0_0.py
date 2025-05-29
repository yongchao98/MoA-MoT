from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 10, 23)

# Subtract 34 days
target_date = start_date - timedelta(days=34)

# Get the full name of the weekday
weekday_name = target_date.strftime('%A')

print(weekday_name)