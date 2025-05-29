# Days in each month from January to March
days_in_january = 31
days_in_february = 28
days_in_march = 31
days_in_april = 5

# Total days from January 1st to April 5th
total_days = days_in_january + days_in_february + days_in_march + days_in_april

# Days of the week starting from Tuesday
days_of_week = ["Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday", "Monday"]

# Calculate the day of the week for April 5th
day_of_week_index = total_days % 7
day_of_week = days_of_week[day_of_week_index]

print(day_of_week)