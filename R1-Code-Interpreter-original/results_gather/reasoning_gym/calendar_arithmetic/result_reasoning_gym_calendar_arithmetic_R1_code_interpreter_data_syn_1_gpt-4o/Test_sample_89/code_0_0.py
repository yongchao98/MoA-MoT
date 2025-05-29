from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 4, 19)

# Number of days to count backward
days_to_subtract = 51

# Calculate the new date
new_date = start_date - timedelta(days=days_to_subtract)

# Get the day of the week
day_of_week = new_date.strftime("%A")

print(day_of_week)