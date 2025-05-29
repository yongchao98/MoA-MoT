from datetime import datetime, timedelta

# Current date
current_date = datetime(2002, 7, 29)

# Add 7 days
new_date = current_date + timedelta(days=7)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")

# Print the formatted date
print(formatted_date)