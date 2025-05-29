# Number of female colleagues in Sarah's circle, including Matilda
female_colleagues_sarah = 5

# Number of female colleagues in John's circle, including Matilda
female_colleagues_john = 4

# Subtract 1 to account for Matilda being counted in both circles
total_female_colleagues = female_colleagues_sarah + female_colleagues_john - 1

# Output the result
print(f"<<<{total_female_colleagues}>>>")