# Days in each month
days_in_january = 31
days_in_february = 28
days_in_march = 28

# Total days from January 1 to March 28
total_days = days_in_january + days_in_february + days_in_march

# Weekdays list starting from Wednesday
weekdays = ["Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday", "Tuesday"]

# Find the weekday for March 28
weekday_index = total_days % 7
weekday_on_march_28 = weekdays[weekday_index]

print(weekday_on_march_28)