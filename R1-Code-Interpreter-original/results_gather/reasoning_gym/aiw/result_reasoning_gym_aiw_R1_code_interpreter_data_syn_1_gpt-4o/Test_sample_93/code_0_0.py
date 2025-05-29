# Number of female colleagues in Linda's circle
female_colleagues_linda = 5

# Number of female colleagues in Charles's circle
female_colleagues_charles = 2

# Matilda is counted in both circles, so we need to find the unique count
# Matilda is one of the female colleagues in both circles
# Therefore, we subtract 1 from each to find the unique count
unique_female_colleagues_linda = female_colleagues_linda - 1
unique_female_colleagues_charles = female_colleagues_charles - 1

# Total unique female colleagues for Matilda
total_female_colleagues = unique_female_colleagues_linda + unique_female_colleagues_charles

print(total_female_colleagues)