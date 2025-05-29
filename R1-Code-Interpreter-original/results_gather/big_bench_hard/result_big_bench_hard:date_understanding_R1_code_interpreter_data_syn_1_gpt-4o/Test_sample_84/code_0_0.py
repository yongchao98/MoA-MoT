from datetime import datetime, timedelta

# Starting date
start_date = datetime.strptime("06/19/2019", "%m/%d/%Y")

# Add 24 hours (1 day)
new_date = start_date + timedelta(days=1)

# Format the new date in MM/DD/YYYY
result = new_date.strftime("%m/%d/%Y")
print(result)