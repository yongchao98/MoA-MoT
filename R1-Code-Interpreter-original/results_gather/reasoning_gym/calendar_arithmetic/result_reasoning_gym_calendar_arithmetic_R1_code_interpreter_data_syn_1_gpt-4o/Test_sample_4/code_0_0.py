import datetime

# January 1, 2023 is a Monday
start_date = datetime.date(2023, 1, 1)
# Calculate the date for January 3, 2023
target_date = start_date + datetime.timedelta(days=2)
# Get the full name of the weekday
weekday_name = target_date.strftime("%A")

print(weekday_name)