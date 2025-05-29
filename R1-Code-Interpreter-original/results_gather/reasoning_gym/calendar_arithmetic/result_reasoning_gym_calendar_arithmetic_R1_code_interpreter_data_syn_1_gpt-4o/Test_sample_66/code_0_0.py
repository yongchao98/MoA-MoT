# Define the weekdays in order
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Given that January 1 is a Friday, find the index of Friday
start_day_index = weekdays.index("Friday")

# Calculate the index of the day 5 days later
days_later = 5
result_day_index = (start_day_index + days_later) % 7

# Find the weekday corresponding to the calculated index
result_day = weekdays[result_day_index]

# Output the result
print(result_day)