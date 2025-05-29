from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 6, 14)

# Number of days to add
days_to_add = 15

# Calculate the new date
new_date = start_date + timedelta(days=days_to_add)

# Get the day of the week
day_of_week = new_date.strftime("%A")

print(day_of_week)