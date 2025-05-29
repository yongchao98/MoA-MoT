from datetime import datetime, timedelta

# Given date
date_str = "07/07/1972"
date_format = "%m/%d/%Y"

# Parse the date
given_date = datetime.strptime(date_str, date_format)

# Subtract 10 days
new_date = given_date - timedelta(days=10)

# Format the new date
result_date = new_date.strftime(date_format)

# Print the result
print(result_date)