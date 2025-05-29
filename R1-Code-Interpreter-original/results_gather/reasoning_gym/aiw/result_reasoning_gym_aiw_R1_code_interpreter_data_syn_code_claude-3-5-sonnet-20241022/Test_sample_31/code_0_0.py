# Create sets of friends
male_friends = set(range(4))  # 4 male friends (0 to 3)
female_friends = set(range(4, 10))  # 6 female friends (4 to 9)

# Let's say William is male friend #0
william_id = 0

# Count William's female friends
williams_female_friends = len(female_friends)

print(f"William has {williams_female_friends} female friends")