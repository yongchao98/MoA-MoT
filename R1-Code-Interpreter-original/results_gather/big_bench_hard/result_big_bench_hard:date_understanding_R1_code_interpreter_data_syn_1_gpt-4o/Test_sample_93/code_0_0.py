from datetime import datetime, timedelta

# Current date
current_date = datetime(1997, 8, 3)

# Subtract 7 days
one_week_ago = current_date - timedelta(days=7)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")

print(formatted_date)