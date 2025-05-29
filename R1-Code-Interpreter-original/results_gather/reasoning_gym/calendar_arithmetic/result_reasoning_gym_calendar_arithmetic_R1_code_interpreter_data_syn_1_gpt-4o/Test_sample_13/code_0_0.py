# Calculate the total number of days from January 1 to March 4
days_in_january = 31 - 1  # Days after January 1
days_in_february = 29  # Leap year
days_in_march = 4  # Up to March 4

total_days = days_in_january + days_in_february + days_in_march

# Calculate the weekday
# Thursday is the 4th day of the week (0: Sunday, 1: Monday, ..., 6: Saturday)
start_weekday = 4  # Thursday
weekday_index = (start_weekday + total_days) % 7

# Map the index to the weekday name
weekdays = ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]
weekday_name = weekdays[weekday_index]

print(weekday_name)