from datetime import datetime, timedelta

# Starting date
start_date = datetime(2020, 3, 20)

# Number of days to add
days_to_add = 176

# Calculate the new date
new_date = start_date + timedelta(days=days_to_add)

# Format the new date in MM/DD/YYYY
formatted_date = new_date.strftime("%m/%d/%Y")

print(formatted_date)