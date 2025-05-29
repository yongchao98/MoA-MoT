from datetime import datetime, timedelta

# Jane quit her job on March 20, 2020
quit_date = datetime(2020, 3, 20)

# Calculate the date 176 days after quitting
days_passed = 176
current_date = quit_date + timedelta(days=days_passed)

# Calculate the date one week ago from the current date
one_week_ago = current_date - timedelta(weeks=1)

# Format the date in MM/DD/YYYY
formatted_date = one_week_ago.strftime("%m/%d/%Y")
print(formatted_date)