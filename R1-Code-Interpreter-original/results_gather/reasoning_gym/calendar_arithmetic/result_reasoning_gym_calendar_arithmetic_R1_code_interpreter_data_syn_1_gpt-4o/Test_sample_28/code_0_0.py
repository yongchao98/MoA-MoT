# Calculate the total number of days from January 1 to March 15
total_days = 31 + 28 + 15

# Calculate the day of the week
# 0 = Wednesday, 1 = Thursday, ..., 6 = Tuesday
day_of_week = (total_days % 7)

# Map the result to the corresponding day of the week
days = ["Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday", "Tuesday"]
result_day = days[day_of_week]

print(result_day)