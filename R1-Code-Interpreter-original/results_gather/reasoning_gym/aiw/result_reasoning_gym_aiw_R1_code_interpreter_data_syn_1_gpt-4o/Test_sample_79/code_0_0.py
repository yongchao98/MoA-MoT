# Number of female colleagues in Linda's circle
female_colleagues_linda = 2

# Number of female colleagues in David's circle
female_colleagues_david = 4

# Matilda is counted in both circles, so we need to find the unique female colleagues
# Matilda is one of the female colleagues in both circles
# Therefore, the unique female colleagues Matilda has are:
unique_female_colleagues = female_colleagues_linda + female_colleagues_david - 1  # Subtracting Matilda herself

print(unique_female_colleagues)