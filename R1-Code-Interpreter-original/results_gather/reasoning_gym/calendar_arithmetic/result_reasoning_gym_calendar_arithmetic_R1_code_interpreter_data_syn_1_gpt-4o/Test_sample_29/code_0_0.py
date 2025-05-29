import calendar

# March 2022
year = 2022
month = 3

# Find all Tuesdays in March 2022
tuesdays = [day for day in range(1, 32) if calendar.weekday(year, month, day) == calendar.TUESDAY]

# Get the last Tuesday
last_tuesday = tuesdays[-1] if tuesdays else -1

print(last_tuesday)