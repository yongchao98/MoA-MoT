# Number of female colleagues in Jennifer's circle (including Matilda)
females_jennifer_circle = 3

# Number of female colleagues in Joseph's circle (including Matilda)
females_joseph_circle = 4

# Matilda is counted in both circles, so we subtract 1 to avoid double-counting her
unique_female_colleagues = females_jennifer_circle + females_joseph_circle - 1

print(unique_female_colleagues)