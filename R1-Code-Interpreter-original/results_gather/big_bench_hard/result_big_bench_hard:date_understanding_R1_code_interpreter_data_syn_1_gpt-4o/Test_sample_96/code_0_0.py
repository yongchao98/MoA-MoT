from datetime import datetime, timedelta

# Date when Jane bought the eggs
start_date = datetime(2017, 5, 9)

# Calculate the date when she ran out of eggs (40 days later)
end_date = start_date + timedelta(days=40)

# Calculate the date 10 days before she ran out of eggs
date_10_days_ago = end_date - timedelta(days=10)

# Format the date in MM/DD/YYYY
formatted_date = date_10_days_ago.strftime("%m/%d/%Y")
print(formatted_date)