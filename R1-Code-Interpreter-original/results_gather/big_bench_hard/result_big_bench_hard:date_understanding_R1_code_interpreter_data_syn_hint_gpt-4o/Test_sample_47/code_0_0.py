from datetime import datetime, timedelta

# Define the date of the golden wedding anniversary
anniversary_date = datetime(2008, 1, 2)

# Calculate the date 10 days before the anniversary
date_10_days_ago = anniversary_date - timedelta(days=10)

# Format the date in MM/DD/YYYY
formatted_date = date_10_days_ago.strftime("%m/%d/%Y")
print(formatted_date)