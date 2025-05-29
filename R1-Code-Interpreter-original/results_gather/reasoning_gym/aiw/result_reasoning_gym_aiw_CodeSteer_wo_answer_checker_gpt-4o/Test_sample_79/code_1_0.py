# Number of female colleagues in Linda's circle (excluding Matilda)
female_colleagues_linda = 2

# Number of female colleagues in David's circle (excluding Matilda)
female_colleagues_david = 4

# Matilda is counted in both circles, so we subtract 1 to avoid double counting
unique_female_colleagues = female_colleagues_linda + female_colleagues_david - 1

# Output the result
print(f"<<<{unique_female_colleagues}>>>")