from datetime import datetime, timedelta

# Define the starting date
start_date = datetime(1899, 12, 31)

# Add 7 days to the starting date
new_date = start_date + timedelta(days=7)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")

# Print the result
print(formatted_date)