# Define Mary's circle
mary_male_colleagues = 3
mary_female_colleagues = 6

# Define Thomas's circle
thomas_male_colleagues = 1
thomas_female_colleagues = 1

# We know Matilda is female and is in both circles
# For Matilda in Mary's circle:
# She is one of the 6 female colleagues of Mary
# So from Matilda's perspective in Mary's circle:
matilda_female_colleagues_from_mary = mary_female_colleagues - 1  # excluding herself

# For Matilda in Thomas's circle:
# She is the female colleague mentioned in Thomas's count
# So from Matilda's perspective in Thomas's circle:
matilda_female_colleagues_from_thomas = 0  # no other females in Thomas's circle except herself

# Total female colleagues of Matilda
total_female_colleagues = matilda_female_colleagues_from_mary

print(f"Matilda has {total_female_colleagues} female colleagues")