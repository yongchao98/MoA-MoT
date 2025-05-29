import datetime

# Define the date
date = datetime.date(2022, 12, 20)

# Get the weekday as an integer
weekday_index = date.weekday()

# Map the integer to the full weekday name
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
weekday_name = weekdays[weekday_index]

# Print the full weekday name
print(weekday_name)