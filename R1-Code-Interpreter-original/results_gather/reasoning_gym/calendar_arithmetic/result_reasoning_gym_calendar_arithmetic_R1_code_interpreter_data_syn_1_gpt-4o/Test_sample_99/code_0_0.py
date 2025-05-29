# Days in each month from January to April
days_in_january = 31
days_in_february = 28
days_in_march = 31
days_in_april = 9

# Total days from January 1st to April 9th
total_days = days_in_january + days_in_february + days_in_march + days_in_april

# Weekdays list starting with Wednesday as January 1st
weekdays = ["Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday", "Tuesday"]

# Calculate the weekday of April 9th
weekday_index = total_days % 7
weekday_of_april_9 = weekdays[weekday_index]

print(weekday_of_april_9)