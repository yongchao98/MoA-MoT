from datetime import datetime, timedelta

# Current date
current_date = datetime(2003, 9, 7)

# Calculate the date 10 days ago
date_10_days_ago = current_date - timedelta(days=10)

# Format the date in MM/DD/YYYY
formatted_date = date_10_days_ago.strftime("%m/%d/%Y")
print(formatted_date)