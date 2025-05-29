import calendar

# October 2022
year = 2022
month = 10

# Find the last day of the month
last_day = calendar.monthrange(year, month)[1]

# Find the weekday of the last day of the month
last_weekday = calendar.weekday(year, month, last_day)

# Calculate the last Monday
# Monday is 0 in Python's calendar module
days_to_last_monday = (last_weekday - 0) % 7
last_monday = last_day - days_to_last_monday

print(last_monday)