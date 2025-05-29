# Linda's circle (including Linda)
linda_males = 3
linda_females = 5  # This includes Matilda

# Charles's circle (including Charles)
charles_males = 6
charles_females = 2  # This includes Matilda

# Since Matilda is in both circles and is female
# Let's calculate Matilda's colleagues:

# In Linda's circle:
# - All males in Linda's circle (3)
# - All females in Linda's circle except herself (5-1 = 4)
matilda_from_linda = linda_males + (linda_females - 1)

# In Charles's circle:
# - All males in Charles's circle (6)
# - All females in Charles's circle except herself (2-1 = 1)
matilda_from_charles = charles_males + (charles_females - 1)

# Total female colleagues of Matilda:
# From Linda's circle: linda_females - 1 (excluding herself)
# From Charles's circle: charles_females - 1 (excluding herself)
matilda_female_colleagues = (linda_females - 1) + (charles_females - 1)

print(f"Matilda has {matilda_female_colleagues} female colleagues")