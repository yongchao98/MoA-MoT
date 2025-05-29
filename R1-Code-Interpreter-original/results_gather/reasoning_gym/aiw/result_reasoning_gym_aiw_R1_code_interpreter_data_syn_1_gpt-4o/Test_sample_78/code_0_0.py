# Number of female colleagues in Susan's circle
female_susan = 4

# Number of female colleagues in Charles's circle
female_charles = 6

# Matilda is counted in both circles, so we subtract 1 to avoid double counting
unique_female_colleagues = female_susan + female_charles - 1

print(unique_female_colleagues)