from datetime import datetime, timedelta

# Current date
current_date = datetime(2004, 5, 4)

# Subtract one year
one_year_ago = current_date.replace(year=current_date.year - 1)

# Format the date in MM/DD/YYYY
formatted_date = one_year_ago.strftime("%m/%d/%Y")
print(formatted_date)