from datetime import datetime, timedelta

# First Monday of 2019
first_day_2019 = datetime(2019, 1, 1)
days_until_monday = (7 - first_day_2019.weekday() + 0) % 7  # 0 is Monday
first_monday_2019 = first_day_2019 + timedelta(days=days_until_monday)

# Date one year before the first Monday of 2019
one_year_ago = first_monday_2019 - timedelta(days=365)

# Format the date in MM/DD/YYYY
formatted_date = one_year_ago.strftime("%m/%d/%Y")
print(formatted_date)