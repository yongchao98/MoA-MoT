from datetime import datetime, timedelta

# Given date
current_date = datetime(2003, 9, 7)

# Subtract 10 days
new_date = current_date - timedelta(days=10)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")
print(formatted_date)