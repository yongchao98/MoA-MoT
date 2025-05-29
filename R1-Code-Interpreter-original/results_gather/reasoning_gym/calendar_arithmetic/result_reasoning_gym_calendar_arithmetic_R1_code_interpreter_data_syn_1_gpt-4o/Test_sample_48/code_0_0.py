# Calculate the total number of days from January 1st to March 11th
days_in_january = 31 - 1  # Days after January 1st
days_in_february = 28
days_in_march = 11

total_days = days_in_january + days_in_february + days_in_march

# Calculate the weekday
# Thursday is the 4th day of the week (0-indexed: 0=Monday, 1=Tuesday, ..., 6=Sunday)
start_weekday = 3  # Thursday

# Find the weekday of March 11
weekday_of_march_11 = (start_weekday + total_days) % 7

# Map the result to the weekday name
weekdays = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
result = weekdays[weekday_of_march_11]

print(result)