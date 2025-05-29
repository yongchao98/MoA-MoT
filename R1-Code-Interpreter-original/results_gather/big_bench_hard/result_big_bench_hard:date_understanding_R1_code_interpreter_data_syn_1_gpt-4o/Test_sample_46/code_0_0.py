from datetime import datetime, timedelta

# Jane's birth date
birth_date = datetime(2001, 2, 28)

# Calculate her 16th birthday
sixteenth_birthday = birth_date.replace(year=birth_date.year + 16)

# Calculate the date 24 hours later
date_24_hours_later = sixteenth_birthday + timedelta(days=1)

# Format the date in MM/DD/YYYY
formatted_date = date_24_hours_later.strftime("%m/%d/%Y")
print(formatted_date)