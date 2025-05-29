from datetime import datetime, timedelta

# Define today's date
today_date = datetime(2011, 1, 22)

# Calculate the date one week ago
one_week_ago = today_date - timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")

print(formatted_date)