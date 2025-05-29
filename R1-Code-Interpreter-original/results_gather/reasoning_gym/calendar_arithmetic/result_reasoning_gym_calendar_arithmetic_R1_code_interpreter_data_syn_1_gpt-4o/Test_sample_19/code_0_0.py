from datetime import datetime

# Create a date object for September 18, 2022
date = datetime(2022, 9, 18)

# Get the full name of the weekday
weekday_name = date.strftime("%A")

# Print the result
print(weekday_name)