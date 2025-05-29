from datetime import datetime, timedelta

# Given date
given_date = datetime(1909, 9, 9)

# Calculate one year ago
one_year_ago = given_date.replace(year=given_date.year - 1)

# Format the date in MM/DD/YYYY
formatted_date = one_year_ago.strftime("%m/%d/%Y")
print(formatted_date)