import calendar

# November 2001
year = 2001
month = 11

# Find the fourth Thursday
month_calendar = calendar.monthcalendar(year, month)
# Thursday is the 4th day of the week (0=Monday, ..., 6=Sunday)
thursday_count = 0
thanksgiving_date = None

for week in month_calendar:
    if week[calendar.THURSDAY] != 0:
        thursday_count += 1
        if thursday_count == 4:
            thanksgiving_date = week[calendar.THURSDAY]
            break

# Calculate the date one month ago
from datetime import date, timedelta

thanksgiving_day = date(year, month, thanksgiving_date)
one_month_ago = thanksgiving_day - timedelta(days=30)

# Output the date in MM/DD/YYYY format
print(one_month_ago.strftime("%m/%d/%Y"))