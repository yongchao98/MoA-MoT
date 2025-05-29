from datetime import datetime, timedelta

# Given date
today_date = datetime(2002, 7, 29)

# Calculate the date one week later
one_week_later = today_date + timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_later.strftime("%m/%d/%Y")
print(formatted_date)