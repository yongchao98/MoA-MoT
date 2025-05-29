# Calculate the total number of days from January 1st to February 25th
total_days = 31 + 25

# Calculate the weekday by finding the remainder when total_days is divided by 7
# Monday is represented as 0, so we add the remainder to 0
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
weekday_index = total_days % 7

# Get the weekday name
weekday_name = weekdays[weekday_index]
print(weekday_name)