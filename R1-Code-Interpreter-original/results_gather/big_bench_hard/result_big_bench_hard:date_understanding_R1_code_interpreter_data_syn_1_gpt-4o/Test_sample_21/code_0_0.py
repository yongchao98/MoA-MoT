from datetime import datetime, timedelta

# Given date
given_date = datetime(2011, 1, 21)

# Subtract 10 days
new_date = given_date - timedelta(days=10)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")

print(formatted_date)