# Number of female colleagues in Linda's circle
female_colleagues_linda = 3

# Number of female colleagues in Richard's circle
female_colleagues_richard = 4

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting her
total_female_colleagues_matilda = female_colleagues_linda + female_colleagues_richard - 1

print(total_female_colleagues_matilda)