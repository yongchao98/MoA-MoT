year = 2048

is_leap_year = (year % 4 == 0) and (year % 100 != 0 or year % 400 == 0)
print("Yes" if is_leap_year else "No")