# Starting day of the week
start_day = "Friday"

# Days of the week in order
days_of_week = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Calculate the index of the starting day
start_index = days_of_week.index(start_day)

# Number of days from January 1 to January 23
days_passed = 22

# Calculate the index of the day of the week for January 23
day_index = (start_index + days_passed) % 7

# Determine the day of the week
day_of_week = days_of_week[day_index]

print(day_of_week)