# Number of female colleagues in Susan's circle
female_colleagues_susan = 4

# Number of female colleagues in Charles's circle, including Matilda
female_colleagues_charles = 6

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting
# Matilda's female colleagues are those in Susan's circle plus those in Charles's circle, excluding Matilda herself
female_colleagues_matilda = female_colleagues_susan + (female_colleagues_charles - 1)

# Output the result
print(f"<<<{female_colleagues_matilda}>>>")