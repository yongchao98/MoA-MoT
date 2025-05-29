from datetime import datetime, timedelta

# Define the last day of the first quarter of 2008
last_day_first_quarter = datetime(2008, 3, 31)

# Calculate the date one week from this day
one_week_later = last_day_first_quarter + timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_later.strftime("%m/%d/%Y")
print(formatted_date)