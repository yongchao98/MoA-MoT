# Calculate the total number of days from January 1 to April 10
total_days = 31 + 28 + 31 + 10

# Calculate the offset from Tuesday
offset = total_days % 7

# Days of the week starting from Tuesday
days_of_week = ["Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday"]

# Determine the day of the week for April 10
day_of_week_april_10 = days_of_week[offset]

print(day_of_week_april_10)