# Calculate the day of the week for February 19
days_from_jan1_to_feb19 = 31 + 19
first_day_of_year = 3  # Thursday is the 4th day of the week if we start counting from Sunday as 0

# Calculate the day of the week for February 19
day_of_week_feb19 = (first_day_of_year + days_from_jan1_to_feb19 - 1) % 7

# Map the result to the corresponding day of the week
days_of_week = ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]
print(days_of_week[day_of_week_feb19])