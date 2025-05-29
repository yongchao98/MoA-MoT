import calendar

# Determine the last day of February 2022
year = 2022
month = 2
last_day = calendar.monthrange(year, month)[1]

# Find the day of the week for the last day of February 2022
# Monday is 0 and Sunday is 6
last_day_weekday = calendar.weekday(year, month, last_day)

# Calculate the day of the month for the last Thursday
# Thursday is 3
days_to_last_thursday = (last_day_weekday - 3) % 7
last_thursday_day = last_day - days_to_last_thursday

print(last_thursday_day)