from datetime import datetime, timedelta

# Define the initial date
initial_date = datetime(1972, 7, 7)

# Subtract 10 days
new_date = initial_date - timedelta(days=10)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")

print(formatted_date)