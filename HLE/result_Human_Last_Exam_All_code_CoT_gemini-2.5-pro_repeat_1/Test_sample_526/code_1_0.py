# You may need to install the 'lunardate' package first.
# You can do this by running the following command in your terminal or command prompt:
# pip install lunardate

import lunardate
from lunardate import LunarDate

# Define the person's Solar birth date.
birth_year = 1980
birth_month = 10
birth_day = 1

# Convert the Solar birth date to the Lunar calendar to find the target Lunar birthday.
lunar_birth_date = LunarDate.fromSolarDate(birth_year, birth_month, birth_day)
target_lunar_month = lunar_birth_date.month
target_lunar_day = lunar_birth_date.day

# Start searching for a match from the year after birth.
current_year = birth_year + 1

while True:
    # For the current year, find out the Lunar date for the Solar date of October 1st.
    current_lunar_date = LunarDate.fromSolarDate(current_year, birth_month, birth_day)

    # Check if the lunar month and day match the person's lunar birthday.
    if (current_lunar_date.month == target_lunar_month and
        current_lunar_date.day == target_lunar_day):
        
        # If they match, we have found the year. Print the final details.
        print(f"The person's Lunar birthday (MM-DD) is: {target_lunar_month}-{target_lunar_day}")
        print(f"In the year {current_year}, the Solar date {current_year}-{birth_month}-{birth_day} is equivalent to the Lunar date {current_lunar_date.month}-{current_lunar_date.day}.")
        print(f"Since {current_lunar_date.month}-{current_lunar_date.day} (Lunar) matches {target_lunar_month}-{target_lunar_day} (Target Lunar), the next year this occurs is:")
        print(current_year)
        
        break  # Exit the loop as we've found the answer.

    # If it's not a match, move to the next year.
    current_year += 1