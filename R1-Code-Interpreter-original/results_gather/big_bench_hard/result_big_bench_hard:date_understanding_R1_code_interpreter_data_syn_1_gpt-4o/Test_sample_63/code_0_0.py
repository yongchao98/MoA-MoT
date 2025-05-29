from datetime import datetime, timedelta

# Original date
original_date = datetime(2002, 7, 29)

# Add 24 hours (1 day)
new_date = original_date + timedelta(days=1)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")
print(formatted_date)