from datetime import datetime, timedelta

# Current date
current_date = datetime(1973, 3, 5)

# Calculate the date one week ago
one_week_ago = current_date - timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")
print(formatted_date)