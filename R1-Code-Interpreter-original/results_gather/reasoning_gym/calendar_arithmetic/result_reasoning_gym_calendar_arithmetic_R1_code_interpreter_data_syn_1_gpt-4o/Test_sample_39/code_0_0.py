# Calculate the total number of days from January 1 to March 21
total_days = 31 + 28 + 21

# Calculate the day of the week
# 0 = Friday, 1 = Saturday, ..., 6 = Thursday
day_of_week = (total_days % 7)

# Map the result to the corresponding day name
days = ["Friday", "Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday"]
result_day = days[day_of_week]

print(result_day)