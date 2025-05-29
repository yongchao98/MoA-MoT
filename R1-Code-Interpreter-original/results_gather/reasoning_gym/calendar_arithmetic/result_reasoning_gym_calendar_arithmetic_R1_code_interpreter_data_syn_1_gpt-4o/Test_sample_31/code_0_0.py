import calendar

# August 2022
year = 2022
month = 8

# Create a calendar for the month
cal = calendar.monthcalendar(year, month)

# Find all Wednesdays in the month
wednesdays = [week[calendar.WEDNESDAY] for week in cal if week[calendar.WEDNESDAY] != 0]

# Get the third Wednesday
third_wednesday = wednesdays[2] if len(wednesdays) >= 3 else -1

print(third_wednesday)