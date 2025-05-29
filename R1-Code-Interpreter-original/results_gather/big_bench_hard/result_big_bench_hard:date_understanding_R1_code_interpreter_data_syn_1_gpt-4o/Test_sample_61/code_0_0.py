from datetime import datetime, timedelta

# Given date for the day before yesterday
day_before_yesterday = datetime.strptime("11/23/1933", "%m/%d/%Y")

# Calculate yesterday's date
yesterday = day_before_yesterday + timedelta(days=1)

# Calculate the date 24 hours later
date_24_hours_later = yesterday + timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = date_24_hours_later.strftime("%m/%d/%Y")
print(formatted_date)