# Sarah's circle
sarah_males = 4
sarah_females = 6  # includes Matilda

# James's circle
james_males = 5
james_females = 3  # includes Matilda

# Since Matilda is counted in both circles' female counts
# and she herself is female, we need to:
# 1. Take all females from Sarah's circle except Matilda (sarah_females - 1)
# 2. Take all females from James's circle except Matilda (james_females - 1)
# These are all Matilda's female colleagues

matilda_female_colleagues = (sarah_females - 1) + (james_females - 1)

print(f"Matilda has {matilda_female_colleagues} female colleagues")