from datetime import datetime, timedelta

# Starting date
start_date = datetime(2022, 10, 24)

# Number of days to add
days_to_add = 24

# Calculate the new date
new_date = start_date + timedelta(days=days_to_add)

# Get the name of the day
day_of_week = new_date.strftime("%A")

print(day_of_week)