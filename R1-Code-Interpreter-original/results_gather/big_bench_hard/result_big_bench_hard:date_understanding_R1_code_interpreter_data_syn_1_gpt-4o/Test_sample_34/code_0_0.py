from datetime import datetime, timedelta

# First Monday of 2019
first_monday_2019 = datetime(2019, 1, 1) + timedelta(days=(7 - 1))  # 6 days after Tuesday

# Date 24 hours later
date_24_hours_later = first_monday_2019 + timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = date_24_hours_later.strftime("%m/%d/%Y")
print(formatted_date)