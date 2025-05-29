# Create sets of friends
male_friends = set(range(4))  # 4 male friends (0 to 3)
female_friends = set(range(4, 10))  # 6 female friends (4 to 9)

# Let's say Bob is male friend #0
bob_id = 0

# Count Bob's female friends (which are all female friends)
bobs_female_friends = len(female_friends)

print(f"Bob has {bobs_female_friends} female friends")