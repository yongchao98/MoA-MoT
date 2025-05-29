from datetime import datetime, timedelta

# Jane's correct date
jane_date = datetime.strptime("06/18/2019", "%m/%d/%Y")

# Calculate the date of yesterday
yesterday_date = jane_date - timedelta(days=1)

# Format the date in MM/DD/YYYY
yesterday_date_str = yesterday_date.strftime("%m/%d/%Y")

print(yesterday_date_str)