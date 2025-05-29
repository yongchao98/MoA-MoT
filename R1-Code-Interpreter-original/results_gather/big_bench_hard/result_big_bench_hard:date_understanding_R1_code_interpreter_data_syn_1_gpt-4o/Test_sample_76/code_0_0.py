from datetime import datetime, timedelta

# Current date is Jane's 3-year work anniversary
current_date = datetime(2019, 12, 1)

# Calculate the date one month ago
one_month_ago = current_date - timedelta(days=30)

# Format the date in MM/DD/YYYY
formatted_date = one_month_ago.strftime("%m/%d/%Y")
print(formatted_date)