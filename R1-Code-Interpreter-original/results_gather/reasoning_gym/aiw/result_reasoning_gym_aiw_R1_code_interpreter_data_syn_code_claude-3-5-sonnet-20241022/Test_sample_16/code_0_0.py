# Create sets of friends
male_friends = set(range(1, 7))  # 6 male friends (1 to 6)
female_friends = set(range(7, 12))  # 5 female friends (7 to 11)

# Let's say Thomas is friend #1
thomas = 1

# Count Thomas's female friends (which are all female friends)
thomas_female_friends = len(female_friends)

print(f"Thomas has {thomas_female_friends} female friends")