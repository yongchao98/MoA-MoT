from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

# Given date
date_str = "12/31/1929"
date_format = "%m/%d/%Y"
given_date = datetime.strptime(date_str, date_format)

# Calculate one month ago
one_month_ago = given_date - relativedelta(months=1)

# Format the result
result_date = one_month_ago.strftime(date_format)
print(result_date)