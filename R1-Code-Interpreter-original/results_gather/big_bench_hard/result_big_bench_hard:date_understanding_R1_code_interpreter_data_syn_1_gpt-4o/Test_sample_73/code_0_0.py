from datetime import datetime, timedelta

# Starting date
start_date = datetime(2012, 1, 31)

# Add 7 days
new_date = start_date + timedelta(days=7)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")
print(formatted_date)