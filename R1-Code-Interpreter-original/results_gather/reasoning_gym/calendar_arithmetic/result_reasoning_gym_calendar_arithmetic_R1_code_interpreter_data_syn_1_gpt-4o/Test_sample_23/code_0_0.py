# Number of days to move forward
days_to_move = 95

# Days of the week starting from Tuesday
days_of_week = ["Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday"]

# Calculate the index of the new day
new_day_index = (days_of_week.index("Tuesday") + days_to_move) % 7

# Get the new day
new_day = days_of_week[new_day_index]

print(new_day)