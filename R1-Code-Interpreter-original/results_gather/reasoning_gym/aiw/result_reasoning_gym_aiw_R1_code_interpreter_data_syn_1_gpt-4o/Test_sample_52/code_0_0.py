# Number of female colleagues in Mary's circle (including Matilda)
females_mary_circle = 5

# Number of female colleagues in Joseph's circle (including Matilda)
females_joseph_circle = 3

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting her
unique_female_colleagues = females_mary_circle + females_joseph_circle - 1

print(unique_female_colleagues)