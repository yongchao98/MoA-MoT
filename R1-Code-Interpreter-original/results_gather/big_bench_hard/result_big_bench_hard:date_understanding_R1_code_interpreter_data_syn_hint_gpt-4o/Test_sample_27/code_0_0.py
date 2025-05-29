from datetime import datetime, timedelta

# Jane quit her job on March 20, 2020
start_date = datetime(2020, 3, 20)

# 176 days have passed since then
current_date = start_date + timedelta(days=176)

# Calculate the date one week ago from the current date
one_week_ago = current_date - timedelta(weeks=1)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")
print(formatted_date)