import datetime

# Define the date
date = datetime.date(2022, 2, 21)

# Get the weekday as an integer
weekday_number = date.weekday()

# Map the integer to the weekday name
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
weekday_name = weekdays[weekday_number]

# Print the result
print(weekday_name)