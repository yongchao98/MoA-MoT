from datetime import datetime, timedelta

# Define the last day of 1899
last_day_1899 = datetime(1899, 12, 31)

# Add 7 days to this date
one_week_later = last_day_1899 + timedelta(days=7)

# Format the result in MM/DD/YYYY
result = one_week_later.strftime("%m/%d/%Y")
print(result)