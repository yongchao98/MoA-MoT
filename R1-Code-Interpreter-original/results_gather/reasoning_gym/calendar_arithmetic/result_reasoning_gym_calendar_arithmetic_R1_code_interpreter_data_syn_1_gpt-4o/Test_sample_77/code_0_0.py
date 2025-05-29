# Given data
start_day = "Tuesday"
days_in_january = 31
days_in_february = 28

# Total days from January 1 to February 28
total_days = days_in_january + days_in_february

# Days of the week
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Find the index of the start day
start_index = weekdays.index(start_day)

# Calculate the index of the day on February 28
feb_28_index = (start_index + total_days - 1) % 7

# Determine the day of the week on February 28
feb_28_day = weekdays[feb_28_index]

print(feb_28_day)