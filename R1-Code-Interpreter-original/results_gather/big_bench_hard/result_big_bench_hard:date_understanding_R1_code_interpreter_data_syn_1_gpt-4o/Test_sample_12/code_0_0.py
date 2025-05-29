from datetime import datetime, timedelta

# Step 1: Calculate the date 176 days after March 20, 2020
start_date = datetime(2020, 3, 20)
days_passed = 176
current_date = start_date + timedelta(days=days_passed)

# Step 2: Determine today's date (which is current_date in this context)
# Step 3: Find the date one year ago from today
one_year_ago = current_date - timedelta(days=365)

# Format the date in MM/DD/YYYY
formatted_date = one_year_ago.strftime("%m/%d/%Y")
print(formatted_date)