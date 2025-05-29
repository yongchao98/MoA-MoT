from datetime import datetime, timedelta

# Given date in UK format (DD/MM/YYYY)
uk_date_str = "02/01/1987"

# Convert to datetime object
uk_date = datetime.strptime(uk_date_str, "%d/%m/%Y")

# Subtract one year
one_year_ago = uk_date.replace(year=uk_date.year - 1)

# Convert to MM/DD/YYYY format
us_date_str = one_year_ago.strftime("%m/%d/%Y")

# Print the result
print(us_date_str)