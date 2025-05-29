from datetime import datetime, timedelta

# Actual date
actual_date = datetime(2002, 3, 12)

# Calculate 10 days ago
date_10_days_ago = actual_date - timedelta(days=10)

# Format the date in MM/DD/YYYY
formatted_date = date_10_days_ago.strftime("%m/%d/%Y")

print(formatted_date)