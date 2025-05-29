from datetime import datetime, timedelta

# Given that yesterday was April 30, 2021
yesterday = datetime(2021, 4, 30)

# Today is one day after yesterday
today = yesterday + timedelta(days=1)

# One week from today
one_week_from_today = today + timedelta(weeks=1)

# Format the date in MM/DD/YYYY
formatted_date = one_week_from_today.strftime("%m/%d/%Y")
print(formatted_date)