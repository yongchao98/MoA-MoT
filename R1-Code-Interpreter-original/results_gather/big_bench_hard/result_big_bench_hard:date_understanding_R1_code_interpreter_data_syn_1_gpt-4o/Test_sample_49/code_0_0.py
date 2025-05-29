from datetime import datetime, timedelta

# Jane and John's wedding date
wedding_date = datetime(1958, 1, 2)

# Calculate the golden anniversary date
golden_anniversary = wedding_date.replace(year=wedding_date.year + 50)

# Calculate the date yesterday
date_yesterday = golden_anniversary - timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date_yesterday = date_yesterday.strftime("%m/%d/%Y")
print(formatted_date_yesterday)