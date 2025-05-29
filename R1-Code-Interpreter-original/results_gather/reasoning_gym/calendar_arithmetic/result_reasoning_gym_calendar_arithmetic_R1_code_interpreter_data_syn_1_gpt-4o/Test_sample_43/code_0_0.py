import calendar

# October 2022
year = 2022
month = 10

# Find the first day of the month and the number of days in the month
first_day_of_month, num_days = calendar.monthrange(year, month)

# Calculate the first Thursday
first_thursday = (3 - first_day_of_month + 7) % 7 + 1

# Calculate the fourth Thursday
fourth_thursday = first_thursday + 3 * 7

# Print the result
print(fourth_thursday)