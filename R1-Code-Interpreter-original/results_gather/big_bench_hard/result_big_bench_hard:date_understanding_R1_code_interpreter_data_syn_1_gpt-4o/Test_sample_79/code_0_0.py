from datetime import datetime, timedelta

# Given date
given_date = datetime(2002, 7, 29)

# Subtract one day
yesterday_date = given_date - timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = yesterday_date.strftime("%m/%d/%Y")
print(formatted_date)