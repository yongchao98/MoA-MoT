from datetime import datetime, timedelta

# Start date when Jane bought the eggs
start_date = datetime(2017, 5, 9)

# Number of days Jane ate eggs
days_eating = 40

# Calculate the date when she ran out of eggs
end_date = start_date + timedelta(days=days_eating)

# Calculate the date of the day before she ran out of eggs
day_before_end_date = end_date - timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = day_before_end_date.strftime("%m/%d/%Y")

print(formatted_date)