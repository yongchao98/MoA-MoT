from datetime import datetime, timedelta

# Start date
start_date = datetime(2020, 3, 20)

# Add 176 days
end_date = start_date + timedelta(days=176)

# Calculate the date of the next day
next_day = end_date + timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = next_day.strftime("%m/%d/%Y")

print(formatted_date)