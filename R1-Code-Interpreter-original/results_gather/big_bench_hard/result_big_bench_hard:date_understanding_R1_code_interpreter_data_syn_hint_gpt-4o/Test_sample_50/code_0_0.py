import datetime

# Date a month ago
date_a_month_ago = datetime.date(1972, 6, 7)

# Check the day of the week
day_of_week = date_a_month_ago.strftime("%A")

print(day_of_week)