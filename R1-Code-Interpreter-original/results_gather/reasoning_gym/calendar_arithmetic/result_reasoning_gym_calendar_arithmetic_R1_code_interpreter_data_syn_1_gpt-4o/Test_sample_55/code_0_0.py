import datetime

# Create a date object for December 7, 2022
date = datetime.date(2022, 12, 7)

# Get the full name of the day of the week
day_of_week = date.strftime("%A")

# Print the result
print(day_of_week)