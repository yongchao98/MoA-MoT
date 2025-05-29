# Sarah's female colleagues
sarah_female_colleagues = 5

# John's female colleagues
john_female_colleagues = 4

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting
unique_female_colleagues = sarah_female_colleagues + john_female_colleagues - 1

print(unique_female_colleagues)