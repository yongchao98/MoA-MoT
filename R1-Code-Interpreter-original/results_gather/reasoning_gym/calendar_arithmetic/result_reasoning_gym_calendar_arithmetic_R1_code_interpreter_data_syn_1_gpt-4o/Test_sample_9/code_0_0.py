import calendar

# Define the year and month
year = 2022
month = 6  # June

# Create a calendar for the month
cal = calendar.monthcalendar(year, month)

# Find the second Tuesday
# The first week may not have a Tuesday, so we need to check the first two weeks
first_week = cal[0]
second_week = cal[1]

# Check if the first Tuesday is in the first week or the second week
if first_week[calendar.TUESDAY] != 0:
    second_tuesday = first_week[calendar.TUESDAY] + 7
else:
    second_tuesday = second_week[calendar.TUESDAY]

print(second_tuesday)