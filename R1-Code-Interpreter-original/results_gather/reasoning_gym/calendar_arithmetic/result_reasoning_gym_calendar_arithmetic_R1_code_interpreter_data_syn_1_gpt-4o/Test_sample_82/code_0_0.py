import datetime

# Define the date
date = datetime.date(2022, 10, 4)

# Get the weekday as an integer
weekday_index = date.weekday()

# Map the integer to the weekday name
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
weekday_name = weekdays[weekday_index]

# Print the result
print(weekday_name)