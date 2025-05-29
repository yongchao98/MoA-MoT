# Create sets of friends
male_friends = set(range(1, 5))  # 4 male friends
female_friends = set(range(5, 10))  # 5 female friends

# Let's say Thomas is male friend #1
thomas = 1

# Count Thomas's female friends
thomas_female_friends = len(female_friends)

print(f"Thomas has {thomas_female_friends} female friends")