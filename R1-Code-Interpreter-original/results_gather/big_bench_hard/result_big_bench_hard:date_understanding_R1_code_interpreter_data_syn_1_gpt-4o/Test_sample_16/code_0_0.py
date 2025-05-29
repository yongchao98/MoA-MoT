from datetime import datetime, timedelta

# Original date
original_date = datetime(1992, 5, 6)

# Add 10 years
current_date = original_date.replace(year=original_date.year + 10)

# Format the date in MM/DD/YYYY
formatted_date = current_date.strftime("%m/%d/%Y")
print(formatted_date)