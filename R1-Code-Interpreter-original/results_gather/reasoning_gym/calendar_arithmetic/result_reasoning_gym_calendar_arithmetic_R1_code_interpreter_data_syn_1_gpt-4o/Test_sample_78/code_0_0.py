import calendar

# October 2022
year = 2022
month = 10

# Find the first day of the month and the number of days in the month
first_weekday, num_days = calendar.monthrange(year, month)

# Calculate the first Thursday
first_thursday = (3 - first_weekday + 7) % 7 + 1

# Calculate the second Thursday
second_thursday = first_thursday + 7

print(second_thursday)