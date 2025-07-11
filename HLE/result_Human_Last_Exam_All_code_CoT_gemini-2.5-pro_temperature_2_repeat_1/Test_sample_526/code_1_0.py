import sys

# We use the 'lunardate' library for calendar conversions.
# This block checks if the library is installed and provides instructions if it's not.
try:
    from lunardate import Solar, Lunar
except ImportError:
    print("Error: The 'lunardate' library is required.")
    print("Please install it by running the following command in your terminal:")
    print("pip install lunardate")
    sys.exit(1)

# Define the person's solar birth date details.
birth_solar_year = 1980
birth_solar_month = 10
birth_solar_day = 1

# Step 1: Find the corresponding lunar birthday.
# Create a Solar date object for the birth date.
solar_birth_date = Solar(birth_solar_year, birth_solar_month, birth_solar_day)
# Convert it to a Lunar date.
lunar_birth_date = solar_birth_date.to_lunar()
lunar_birth_month = lunar_birth_date.month
lunar_birth_day = lunar_birth_date.day

print(f"The solar birth date is {birth_solar_year}-{birth_solar_month}-{birth_solar_day}.")
print(f"This corresponds to the lunar birth date of month {lunar_birth_month}, day {lunar_birth_day}.")
print("-" * 20)
print(f"Now, we will find the next year 'YYYY' where the lunar date {lunar_birth_month}-{lunar_birth_day} falls on the solar date {birth_solar_month}-{birth_solar_day}.")

# Step 2: Iterate through the years starting from the year after birth.
found_year = None
for age in range(1, 100): # Loop through ages 1 to 99.
    # The lunar year for the birthday anniversary.
    current_lunar_year = birth_solar_year + age
    
    try:
        # Step 3: Create the lunar anniversary date for the current year.
        anniversary_lunar_date = Lunar(current_lunar_year, lunar_birth_month, lunar_birth_day)
        
        # Convert it back to a solar date.
        anniversary_solar_date = anniversary_lunar_date.to_solar()
        
        # Step 4: Check if the solar date matches the original solar birthday's month and day.
        if (anniversary_solar_date.month == birth_solar_month and
            anniversary_solar_date.day == birth_solar_day):
            
            # If it matches, we have found the year.
            found_year = anniversary_solar_date.year
            break # Exit the loop since we found the first occurrence.
            
    except ValueError:
        # This handles cases where a lunar date is invalid in a given year
        # (e.g., day 30 in a month with only 29 days).
        continue

# Step 5: Print the final result.
if found_year:
    print(f"\nThe next year this occurs is: {found_year}")
else:
    print("\nA matching year was not found in the next 100 years.")
