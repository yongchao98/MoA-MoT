# Calculate the total number of days from January 1 to March 6
days_in_january = 30  # From January 2 to January 31
days_in_february = 28
days_in_march = 6

total_days = days_in_january + days_in_february + days_in_march

# Calculate the weekday
# January 1 is a Monday, so we start counting from 0 (Monday)
weekday_index = total_days % 7

# List of weekdays starting from Monday
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Determine the weekday
weekday_on_march_6 = weekdays[weekday_index]

print(weekday_on_march_6)