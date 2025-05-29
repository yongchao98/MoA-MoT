# Given that January 1 is a Monday, calculate the weekday for January 3
# January 1 is a Monday, so January 3 should be two days after Monday

# Days of the week starting from Monday
days_of_week = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# January 1 is a Monday, which is index 0
start_day_index = 0

# Calculate the index for January 3
target_day_index = (start_day_index + 2) % 7

# Get the weekday name
weekday_name = days_of_week[target_day_index]

print(weekday_name)