from datetime import datetime, timedelta

# Today's date is their 5-year anniversary
today_date = datetime(1963, 1, 2)

# Calculate the date one week from today
one_week_later = today_date + timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_later.strftime("%m/%d/%Y")
print(formatted_date)