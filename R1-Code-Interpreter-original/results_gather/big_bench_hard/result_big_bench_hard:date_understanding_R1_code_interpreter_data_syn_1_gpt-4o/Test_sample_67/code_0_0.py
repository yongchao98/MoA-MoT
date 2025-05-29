from datetime import datetime, timedelta

# Today's date is September 8th, 2021
today_date = datetime(2021, 9, 8)

# Calculate the date one week ago
one_week_ago = today_date - timedelta(weeks=1)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")
print(formatted_date)