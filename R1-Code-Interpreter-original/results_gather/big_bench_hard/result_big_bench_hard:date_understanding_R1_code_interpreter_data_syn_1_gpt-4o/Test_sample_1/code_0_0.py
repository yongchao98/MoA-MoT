from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

# Given date in UK format (DD/MM/YYYY)
uk_date_str = "02/01/1987"
uk_date = datetime.strptime(uk_date_str, "%d/%m/%Y")

# Calculate one month ago
one_month_ago = uk_date - relativedelta(months=1)

# Convert to MM/DD/YYYY format
us_date_str = one_month_ago.strftime("%m/%d/%Y")
print(us_date_str)