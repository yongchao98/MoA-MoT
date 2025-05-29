from datetime import datetime, timedelta

# Define the deadline date
deadline_date = datetime(2021, 6, 1)

# Calculate the date 10 days before the deadline
date_10_days_ago = deadline_date - timedelta(days=10)

# Format the date in MM/DD/YYYY format
formatted_date = date_10_days_ago.strftime("%m/%d/%Y")

# Print the formatted date
print(formatted_date)