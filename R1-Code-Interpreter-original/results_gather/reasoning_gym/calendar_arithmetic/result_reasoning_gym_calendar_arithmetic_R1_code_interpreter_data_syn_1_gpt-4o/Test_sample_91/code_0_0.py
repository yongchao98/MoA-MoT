# Starting day is Tuesday
starting_day = "Tuesday"

# Days of the week in order
days_of_week = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

# Find the index of the starting day
start_index = days_of_week.index(starting_day)

# Number of days to jump
days_to_jump = 99

# Calculate the index of the new day
new_index = (start_index + days_to_jump) % 7

# Find the new day
new_day = days_of_week[new_index]

print(new_day)