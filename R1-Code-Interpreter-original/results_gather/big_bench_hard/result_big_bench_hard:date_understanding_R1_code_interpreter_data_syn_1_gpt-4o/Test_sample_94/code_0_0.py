from datetime import datetime, timedelta

# Given date
date_str = "03/02/1966"
date_format = "%m/%d/%Y"

# Convert string to datetime object
date_obj = datetime.strptime(date_str, date_format)

# Add one day (24 hours)
new_date_obj = date_obj + timedelta(days=1)

# Convert back to string in MM/DD/YYYY format
new_date_str = new_date_obj.strftime(date_format)

print(new_date_str)