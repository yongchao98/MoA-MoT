import calendar

# Define the year and month
year = 2022
month = 5  # May

# Create a calendar for the month
cal = calendar.monthcalendar(year, month)

# Find the third Friday
friday_count = 0
third_friday = -1

for week in cal:
    if week[calendar.FRIDAY] != 0:  # Check if there's a Friday in this week
        friday_count += 1
        if friday_count == 3:
            third_friday = week[calendar.FRIDAY]
            break

print(third_friday)