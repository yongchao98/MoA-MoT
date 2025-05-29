import calendar

# Step 1: Identify the day of the week for April 1, 2022
year = 2022
month = 4
first_day_of_month, days_in_month = calendar.monthrange(year, month)

# Step 2: Find the first Tuesday in April 2022
# 0 = Monday, 1 = Tuesday, ..., 6 = Sunday
first_tuesday = 1 + (1 - first_day_of_month) % 7

# Step 3: Calculate the date of the third Tuesday
third_tuesday = first_tuesday + 14

# Print the result
print(third_tuesday)