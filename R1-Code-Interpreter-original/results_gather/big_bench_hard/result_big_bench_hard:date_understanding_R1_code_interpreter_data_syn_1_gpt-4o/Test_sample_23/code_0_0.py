from datetime import datetime, timedelta

# Today's date based on the information given
today = datetime(2014, 12, 30)

# Calculate the date one week ago
one_week_ago = today - timedelta(weeks=1)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")
print(formatted_date)