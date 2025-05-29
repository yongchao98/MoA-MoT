import calendar

# November 2022
year = 2022
month = 11

# Find all Fridays in November 2022
fridays = [day for day in range(1, calendar.monthrange(year, month)[1] + 1)
           if calendar.weekday(year, month, day) == calendar.FRIDAY]

# The last Friday
last_friday = fridays[-1] if fridays else -1

print(last_friday)