from datetime import datetime

# Define the date
date_str = "2022-04-27"

# Convert the string to a datetime object
date_obj = datetime.strptime(date_str, "%Y-%m-%d")

# Get the weekday as an integer
weekday_int = date_obj.weekday()

# Map the integer to the weekday name
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
weekday_name = weekdays[weekday_int]

# Print the result
print(weekday_name)