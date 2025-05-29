# Define Jessica's circle (excluding Jessica herself)
jessica_males = 6
jessica_females = 3

# Define Charles's circle (excluding Charles himself)
charles_males = 1
charles_females = 1

# We know Matilda is female and is in both circles
# So in Jessica's circle, out of 3 females, one is Matilda
# And in Charles's circle, out of 1 female, it's Matilda

# Matilda's female colleagues would be:
# - All females in Jessica's circle (except herself)
# - All females in Charles's circle (except herself)
# - No duplicates as Matilda is the only connection

# Calculate Matilda's female colleagues
matilda_female_colleagues = (jessica_females - 1) + (charles_females - 1)

print(f"Matilda has {matilda_female_colleagues} female colleagues")