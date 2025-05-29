from datetime import datetime, timedelta

# Given date in UK format
uk_date_str = "02/01/1987"

# Convert to datetime object
uk_date = datetime.strptime(uk_date_str, "%d/%m/%Y")

# Find the previous day
previous_day = uk_date - timedelta(days=1)

# Convert to US format
us_date_str = previous_day.strftime("%m/%d/%Y")

# Print the result
print(us_date_str)