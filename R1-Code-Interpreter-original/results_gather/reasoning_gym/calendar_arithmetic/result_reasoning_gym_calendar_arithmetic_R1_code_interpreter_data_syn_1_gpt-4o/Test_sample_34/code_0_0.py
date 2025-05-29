from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 5, 14)

# Calculate the date after 32 days
future_date = start_date + timedelta(days=32)

# Get the weekday name
weekday_name = future_date.strftime('%A')

print(weekday_name)