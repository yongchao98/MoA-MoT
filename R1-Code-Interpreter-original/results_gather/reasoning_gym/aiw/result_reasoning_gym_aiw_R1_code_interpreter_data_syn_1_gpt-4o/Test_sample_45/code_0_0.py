# Define the number of female colleagues in each circle
barbara_female_colleagues = 2
michael_female_colleagues = 3

# Matilda is counted in both circles, so we subtract 1 to avoid double counting
unique_female_colleagues = barbara_female_colleagues + michael_female_colleagues - 1

# Output the number of female colleagues Matilda has
print(unique_female_colleagues)