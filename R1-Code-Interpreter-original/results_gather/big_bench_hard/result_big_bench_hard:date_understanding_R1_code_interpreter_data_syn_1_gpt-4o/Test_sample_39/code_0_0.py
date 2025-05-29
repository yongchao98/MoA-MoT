from datetime import datetime, timedelta

# Given date
start_date = datetime(1969, 4, 19)

# Add 7 days
new_date = start_date + timedelta(days=7)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")
print(formatted_date)