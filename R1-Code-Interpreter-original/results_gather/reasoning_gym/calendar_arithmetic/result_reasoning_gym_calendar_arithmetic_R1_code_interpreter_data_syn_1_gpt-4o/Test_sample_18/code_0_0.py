import calendar

# Define the year and month
year = 2022
month = 5

# Create a calendar for May 2022
cal = calendar.monthcalendar(year, month)

# Find the third Saturday
saturdays = [week[calendar.SATURDAY] for week in cal if week[calendar.SATURDAY] != 0]

# Get the third Saturday
third_saturday = saturdays[2] if len(saturdays) >= 3 else -1

print(third_saturday)