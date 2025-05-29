# Number of female colleagues in Susan's circle
female_colleagues_susan = 5

# Number of female colleagues in James's circle
female_colleagues_james = 6

# Since Matilda is counted in both circles, we subtract 1 to avoid double-counting her
total_female_colleagues_matilda = female_colleagues_susan + female_colleagues_james - 1

print(total_female_colleagues_matilda)