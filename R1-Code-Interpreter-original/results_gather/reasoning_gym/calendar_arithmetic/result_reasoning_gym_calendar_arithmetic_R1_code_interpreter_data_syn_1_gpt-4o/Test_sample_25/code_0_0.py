# Given information
current_day = "Sunday"
days_ago = 78

# Days of the week in order
days_of_week = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Find the index of the current day
current_index = days_of_week.index(current_day)

# Calculate the index of the day 78 days ago
days_ago_index = (current_index - days_ago) % 7

# Find the day of the week 78 days ago
day_78_days_ago = days_of_week[days_ago_index]

print(day_78_days_ago)