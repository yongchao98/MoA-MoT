from datetime import datetime, timedelta

# Given date in UK format
uk_date_str = "02/01/1987"

# Convert to datetime object
uk_date = datetime.strptime(uk_date_str, "%d/%m/%Y")

# Calculate the next day
next_day = uk_date + timedelta(days=1)

# Convert to US format
us_date_str = next_day.strftime("%m/%d/%Y")

print(us_date_str)