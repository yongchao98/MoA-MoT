# Number of female colleagues in Jennifer's circle
female_colleagues_jennifer = 3

# Number of female colleagues in Joseph's circle
female_colleagues_joseph = 4

# Matilda is counted in both circles, so we subtract her from Joseph's circle
# to find the unique female colleagues Matilda has
unique_female_colleagues_matilda = female_colleagues_joseph - 1

# Output the final answer
print(f"<<<{unique_female_colleagues_matilda}>>>")