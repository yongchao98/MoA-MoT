from datetime import datetime, timedelta

# Define the date
date_today = datetime(1972, 7, 7)

# Calculate tomorrow's date
date_tomorrow = date_today + timedelta(days=1)

# Print tomorrow's date
print(date_tomorrow.strftime("%m/%d/%Y"))