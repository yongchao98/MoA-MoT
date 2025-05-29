import datetime

# Create a date object for May 13, 2022
date = datetime.date(2022, 5, 13)

# Get the day of the week as an integer (Monday=0, Sunday=6)
day_of_week_index = date.weekday()

# Map the integer to the corresponding weekday name
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
day_of_week_name = weekdays[day_of_week_index]

# Print the result
print(day_of_week_name)