from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 5, 23)

# Number of days to add
days_to_add = 59

# Calculate the future date
future_date = start_date + timedelta(days=days_to_add)

# Get the full weekday name
future_weekday = future_date.strftime("%A")

print(future_weekday)