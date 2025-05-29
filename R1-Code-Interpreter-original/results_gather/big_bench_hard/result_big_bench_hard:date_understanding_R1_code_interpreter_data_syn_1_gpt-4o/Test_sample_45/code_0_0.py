from datetime import datetime, timedelta

# Today's date after the delay
today_date = datetime.strptime("06/02/1943", "%m/%d/%Y")

# Calculate yesterday's date
yesterday_date = today_date - timedelta(days=1)

# Format the date in MM/DD/YYYY
yesterday_date_str = yesterday_date.strftime("%m/%d/%Y")
print(yesterday_date_str)