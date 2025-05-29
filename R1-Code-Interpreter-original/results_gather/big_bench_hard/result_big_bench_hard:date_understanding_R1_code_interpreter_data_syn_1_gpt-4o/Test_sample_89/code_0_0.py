from datetime import datetime, timedelta

# Current date (today) after the delay
current_date = datetime.strptime("06/02/1943", "%m/%d/%Y")

# Calculate tomorrow's date
tomorrow_date = current_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
print(tomorrow_date.strftime("%m/%d/%Y"))