from datetime import datetime, timedelta

# Start date
start_date = datetime(2020, 3, 20)

# Add 176 days
end_date = start_date + timedelta(days=176)

# Calculate the date for tomorrow
tomorrow_date = end_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = tomorrow_date.strftime("%m/%d/%Y")
print(formatted_date)