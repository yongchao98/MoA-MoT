# Number of female colleagues in Patricia's circle
female_colleagues_patricia = 6

# Number of female colleagues in William's circle
female_colleagues_william = 6

# Matilda is counted in both circles, so subtract her once
matilda_female_colleagues = female_colleagues_patricia + female_colleagues_william - 1

# Output the result
print(f"<<<{matilda_female_colleagues}>>>")