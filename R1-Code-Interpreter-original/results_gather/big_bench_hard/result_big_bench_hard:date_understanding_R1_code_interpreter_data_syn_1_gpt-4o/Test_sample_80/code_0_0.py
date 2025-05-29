import calendar

# November 2001
year = 2001
month = 11

# Find the fourth Thursday
month_calendar = calendar.monthcalendar(year, month)
# A month calendar is a matrix where each row represents a week
# and each column represents a day of the week (0=Monday, 6=Sunday)

# Find the first Thursday
first_thursday = None
for week in month_calendar:
    if week[calendar.THURSDAY] != 0:
        first_thursday = week[calendar.THURSDAY]
        break

# Calculate the fourth Thursday
fourth_thursday = first_thursday + 3 * 7

# Calculate the date 24 hours later
from datetime import datetime, timedelta

thanksgiving_date = datetime(year, month, fourth_thursday)
next_day = thanksgiving_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
next_day_formatted = next_day.strftime("%m/%d/%Y")
print(next_day_formatted)