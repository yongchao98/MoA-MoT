# Number of female colleagues in Linda's circle
female_colleagues_linda = 5

# Number of female colleagues in Charles's circle
female_colleagues_charles = 2

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting her
total_female_colleagues = female_colleagues_linda + female_colleagues_charles - 1

# Output the result
print(f"<<<{total_female_colleagues}>>>")