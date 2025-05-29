from datetime import datetime, timedelta

# Define the date
current_date = datetime(1909, 9, 9)

# Subtract one day
previous_date = current_date - timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = previous_date.strftime("%m/%d/%Y")
print(formatted_date)