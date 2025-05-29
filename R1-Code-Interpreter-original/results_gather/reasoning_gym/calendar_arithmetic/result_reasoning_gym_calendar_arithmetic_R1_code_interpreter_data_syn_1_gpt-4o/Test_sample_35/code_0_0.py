import datetime

# Starting date
start_date = datetime.date(2022, 2, 25)

# Number of days to add
days_to_add = 53

# Calculate the new date
new_date = start_date + datetime.timedelta(days=days_to_add)

# Get the weekday name
weekday_name = new_date.strftime("%A")

print(weekday_name)