from datetime import datetime, timedelta

# Define the palindrome date
palindrome_date_str = "02/04/2020"
palindrome_date = datetime.strptime(palindrome_date_str, "%m/%d/%Y")

# Calculate the next day
next_day = palindrome_date + timedelta(days=1)
next_day_str = next_day.strftime("%m/%d/%Y")

print(next_day_str)