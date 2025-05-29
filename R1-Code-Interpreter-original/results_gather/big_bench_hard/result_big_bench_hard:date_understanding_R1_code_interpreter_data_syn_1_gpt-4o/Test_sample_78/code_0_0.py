from datetime import datetime, timedelta

# Starting date
start_date = datetime(2012, 1, 31)

# Add 24 hours
next_day = start_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = next_day.strftime("%m/%d/%Y")
print(formatted_date)