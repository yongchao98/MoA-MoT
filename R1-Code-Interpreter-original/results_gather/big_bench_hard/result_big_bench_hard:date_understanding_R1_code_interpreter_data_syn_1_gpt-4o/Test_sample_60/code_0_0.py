from datetime import datetime, timedelta

# Current date
current_date = datetime(1969, 4, 19)

# Add 24 hours (1 day)
new_date = current_date + timedelta(days=1)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")
print(formatted_date)