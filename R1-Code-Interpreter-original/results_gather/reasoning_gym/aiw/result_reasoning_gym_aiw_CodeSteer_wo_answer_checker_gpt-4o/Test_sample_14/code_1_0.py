# Number of female colleagues in Susan's circle
female_colleagues_susan = 5

# Number of female colleagues in James's circle
female_colleagues_james = 6

# Matilda is counted in both circles, so subtract her once
matilda_female_colleagues = female_colleagues_susan + female_colleagues_james - 1

# Output the result
print(f"<<<{matilda_female_colleagues}>>>")