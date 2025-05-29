# Define the number of female colleagues in each circle
female_colleagues_susan = 4
female_colleagues_charles = 6

# Matilda is part of both circles, so we need to account for her being counted twice
# We subtract 1 from the total count of female colleagues in Charles's circle to avoid double-counting Matilda
# Matilda's female colleagues are those in Susan's circle plus those in Charles's circle, excluding Matilda herself
female_colleagues_matilda = female_colleagues_susan + (female_colleagues_charles - 1)

# Output the result
print(f"<<<{female_colleagues_matilda}>>>")