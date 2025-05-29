# Create sets of friends
male_friends = set(range(1, 5))  # 4 male friends
female_friends = set(range(5, 10))  # 5 female friends

# Robert is one of the male friends (let's say he's number 1)
robert = 1

# Count Robert's female friends
roberts_female_friends = female_friends  # Since all are friends with each other

print(f"Robert has {len(roberts_female_friends)} female friends")