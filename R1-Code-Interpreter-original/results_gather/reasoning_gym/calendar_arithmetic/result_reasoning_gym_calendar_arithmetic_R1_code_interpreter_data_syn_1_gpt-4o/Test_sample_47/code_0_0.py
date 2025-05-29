import datetime

# Create a date object for 08/12/2022
date = datetime.datetime.strptime('08/12/2022', '%m/%d/%Y')

# Get the day of the week as an integer
day_of_week = date.weekday()

# Map the integer to the corresponding weekday name
weekdays = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']
weekday_name = weekdays[day_of_week]

# Print the result
print(weekday_name)