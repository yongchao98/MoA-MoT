# Mary's circle (excluding Mary)
mary_males = 1
mary_females = 5  # includes Matilda

# Joseph's circle (excluding Joseph)
joseph_males = 6
joseph_females = 3  # includes Matilda

# Since Matilda is counted in both circles, we need to:
# 1. Count female colleagues from Mary's circle (excluding Matilda)
# 2. Count female colleagues from Joseph's circle (excluding Matilda)
# 3. Remove duplicates (there are none as Matilda is the only overlap)

# Female colleagues from Mary's circle (excluding Matilda)
matilda_females_from_mary = mary_females - 1  # -1 because we exclude Matilda herself

# Female colleagues from Joseph's circle (excluding Matilda)
matilda_females_from_joseph = joseph_females - 1  # -1 because we exclude Matilda herself

# Total female colleagues of Matilda
total_female_colleagues = matilda_females_from_mary + matilda_females_from_joseph

print(total_female_colleagues)