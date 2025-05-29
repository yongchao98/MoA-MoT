from datetime import datetime, timedelta

# Given date
today_date = datetime(1972, 7, 9)

# Calculate one week ago
one_week_ago = today_date - timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")

print(formatted_date)