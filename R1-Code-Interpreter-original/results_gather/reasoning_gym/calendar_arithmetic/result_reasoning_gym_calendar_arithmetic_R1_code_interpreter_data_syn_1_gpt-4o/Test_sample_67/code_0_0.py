import datetime

# Define the year and month
year = 2022
month = 7

# Find the first day of the month
first_day_of_month = datetime.date(year, month, 1)

# Calculate the first Friday of the month
first_friday = first_day_of_month + datetime.timedelta((4 - first_day_of_month.weekday() + 7) % 7)

# Calculate the second Friday
second_friday = first_friday + datetime.timedelta(days=7)

# Output the day of the month for the second Friday
print(second_friday.day)