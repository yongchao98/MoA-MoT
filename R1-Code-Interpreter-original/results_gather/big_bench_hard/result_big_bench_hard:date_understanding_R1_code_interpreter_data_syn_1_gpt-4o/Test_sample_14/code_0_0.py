from datetime import datetime, timedelta

# Today's date
today_date = datetime.strptime("11/11/2019", "%m/%d/%Y")

# Calculate the date one week from today
one_week_later = today_date + timedelta(days=7)

# Format the result in MM/DD/YYYY
result_date = one_week_later.strftime("%m/%d/%Y")
print(result_date)