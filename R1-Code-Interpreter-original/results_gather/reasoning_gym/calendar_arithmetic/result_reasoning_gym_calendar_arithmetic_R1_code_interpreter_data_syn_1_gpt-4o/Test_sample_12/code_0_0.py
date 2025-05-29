import calendar

# September 2022
year = 2022
month = 9

# Find the last Saturday
last_saturday = -1
for day in range(30, 0, -1):
    if calendar.weekday(year, month, day) == calendar.SATURDAY:
        last_saturday = day
        break

print(last_saturday)