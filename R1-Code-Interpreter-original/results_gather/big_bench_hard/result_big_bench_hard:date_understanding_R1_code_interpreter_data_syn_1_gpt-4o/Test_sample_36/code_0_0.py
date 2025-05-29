from datetime import datetime, timedelta

# Start date
start_date = datetime(2020, 3, 20)

# Add 176 days
date_after_176_days = start_date + timedelta(days=176)

# Add 24 hours (1 day)
final_date = date_after_176_days + timedelta(days=1)

# Format the final date in MM/DD/YYYY
formatted_date = final_date.strftime("%m/%d/%Y")

print(formatted_date)