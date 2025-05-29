from datetime import datetime, timedelta

# Today's date
today_date = datetime.strptime("10/15/1924", "%m/%d/%Y")

# Calculate the date one week from today
one_week_later = today_date + timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_later.strftime("%m/%d/%Y")
print(formatted_date)