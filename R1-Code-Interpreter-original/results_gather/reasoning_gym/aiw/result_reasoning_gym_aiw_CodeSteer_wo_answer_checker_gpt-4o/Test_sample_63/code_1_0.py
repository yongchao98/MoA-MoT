# Number of female colleagues in Sarah's circle
sarah_female_colleagues = 4

# Number of female colleagues in Thomas's circle
thomas_female_colleagues = 2

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting her
unique_female_colleagues = sarah_female_colleagues + thomas_female_colleagues - 1

# Output the result
print(f"<<<{unique_female_colleagues}>>>")