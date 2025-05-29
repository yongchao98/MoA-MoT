year = 1970

# Check if the year is a leap year
is_leap_year = (year % 4 == 0) and (year % 100 != 0 or year % 400 == 0)

# Output the result
print("Yes" if is_leap_year else "No")