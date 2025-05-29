# Number of female colleagues in Jennifer's circle
female_colleagues_jennifer = 6

# Number of female colleagues in William's circle
female_colleagues_william = 3

# Matilda is counted in both circles, so we subtract 1 to avoid double counting
female_colleagues_matilda = female_colleagues_jennifer + female_colleagues_william - 1

# Output the result
print(f"<<<{female_colleagues_matilda}>>>")