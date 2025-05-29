# Number of female colleagues in Mary's circle (including Matilda)
female_colleagues_mary = 3

# Number of female colleagues in Thomas's circle (including Matilda)
female_colleagues_thomas = 1

# Matilda is counted in both circles, so we subtract 1 to exclude her
unique_female_colleagues = (female_colleagues_mary - 1) + (female_colleagues_thomas - 1)

# Since Matilda is counted twice, we add 1 back to account for her being a unique connection
unique_female_colleagues += 1

print(unique_female_colleagues)