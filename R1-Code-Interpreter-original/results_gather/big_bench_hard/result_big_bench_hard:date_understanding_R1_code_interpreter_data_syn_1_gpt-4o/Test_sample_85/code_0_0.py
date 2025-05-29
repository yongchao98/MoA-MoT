from datetime import datetime, timedelta

# Define today's date
today_date = datetime(2011, 1, 22)

# Calculate the date one week from today
one_week_later = today_date + timedelta(weeks=1)

# Format the date in MM/DD/YYYY
formatted_date = one_week_later.strftime("%m/%d/%Y")

# Print the formatted date
print(formatted_date)