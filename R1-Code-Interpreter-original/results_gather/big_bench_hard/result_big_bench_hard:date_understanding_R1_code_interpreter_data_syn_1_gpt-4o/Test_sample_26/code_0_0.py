from datetime import datetime, timedelta

# Today's date based on the given deadline
today = datetime(2021, 5, 30)

# Calculate the date one week ago
one_week_ago = today - timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")
print(formatted_date)