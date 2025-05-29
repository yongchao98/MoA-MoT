from datetime import datetime, timedelta

# Anniversary date
anniversary_date = datetime(2008, 1, 2)

# Calculate the date 10 days before
date_10_days_before = anniversary_date - timedelta(days=10)

# Format the date in MM/DD/YYYY
formatted_date = date_10_days_before.strftime("%m/%d/%Y")
print(formatted_date)