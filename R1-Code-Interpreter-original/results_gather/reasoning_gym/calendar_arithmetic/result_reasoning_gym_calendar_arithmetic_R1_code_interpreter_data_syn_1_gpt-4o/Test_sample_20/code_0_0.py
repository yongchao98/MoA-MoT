# Calculate the total number of days from January 1st to April 3rd
total_days = 31 + 28 + 31 + 3

# Calculate the weekday by finding the remainder when total_days is divided by 7
# Monday is 0, Tuesday is 1, ..., Sunday is 6
weekday_index = (total_days - 1) % 7

# List of weekdays
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Determine the weekday
weekday = weekdays[weekday_index]

print(weekday)