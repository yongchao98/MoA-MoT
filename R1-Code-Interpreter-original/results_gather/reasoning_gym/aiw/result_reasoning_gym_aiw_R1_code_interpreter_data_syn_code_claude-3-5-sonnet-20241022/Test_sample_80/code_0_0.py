# Create sets of friends
male_friends = set(range(1, 6))  # 5 male friends (1 to 5), let's say Charles is 1
female_friends = set(range(6, 9))  # 3 female friends (6 to 8)

# Charles's ID is 1
charles_id = 1

# Get Charles's female friends (intersection with female_friends set)
charles_female_friends = female_friends

print(f"Number of Charles's female friends: {len(charles_female_friends)}")