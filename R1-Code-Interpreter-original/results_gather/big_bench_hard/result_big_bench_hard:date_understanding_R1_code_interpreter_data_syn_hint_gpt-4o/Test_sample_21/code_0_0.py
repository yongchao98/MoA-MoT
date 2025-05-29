from datetime import datetime, timedelta

# Define the date
date_str = "01/21/2011"
date_format = "%m/%d/%Y"

# Convert string to datetime object
date_obj = datetime.strptime(date_str, date_format)

# Subtract 10 days
new_date_obj = date_obj - timedelta(days=10)

# Convert back to string in the desired format
new_date_str = new_date_obj.strftime(date_format)

print(new_date_str)