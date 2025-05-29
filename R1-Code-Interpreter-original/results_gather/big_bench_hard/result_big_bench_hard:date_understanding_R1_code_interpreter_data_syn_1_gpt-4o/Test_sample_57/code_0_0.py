from datetime import datetime, timedelta

# Jane's 16th birthday
birthday_16 = datetime(2016, 2, 29)

# Calculate one year ago
one_year_ago = birthday_16 - timedelta(days=365)

# Output the date in MM/DD/YYYY format
print(one_year_ago.strftime("%m/%d/%Y"))