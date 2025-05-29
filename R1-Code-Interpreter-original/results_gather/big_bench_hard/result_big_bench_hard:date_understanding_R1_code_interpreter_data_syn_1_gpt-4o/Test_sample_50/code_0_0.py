from datetime import datetime, timedelta

# Given date
given_date = datetime(1972, 7, 9)

# Calculate the date one month ago
one_month_ago = given_date - timedelta(days=30)

# Print the result in MM/DD/YYYY format
print(one_month_ago.strftime("%m/%d/%Y"))