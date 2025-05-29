# Number of female colleagues in Margaret's circle
female_colleagues_margaret = 3

# Number of female colleagues in Thomas's circle
female_colleagues_thomas = 6

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting her
unique_female_colleagues = female_colleagues_margaret + female_colleagues_thomas - 1

# Output the result
print(f"<<<{unique_female_colleagues}>>>")