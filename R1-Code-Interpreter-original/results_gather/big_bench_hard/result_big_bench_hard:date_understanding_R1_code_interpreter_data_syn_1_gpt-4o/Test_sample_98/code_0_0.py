from datetime import datetime, timedelta

# Jane's 16th birthday
birthday_16 = datetime(2017, 2, 28)

# Calculate the date 10 days before
date_10_days_before = birthday_16 - timedelta(days=10)

# Format the date in MM/DD/YYYY
formatted_date = date_10_days_before.strftime("%m/%d/%Y")
print(formatted_date)