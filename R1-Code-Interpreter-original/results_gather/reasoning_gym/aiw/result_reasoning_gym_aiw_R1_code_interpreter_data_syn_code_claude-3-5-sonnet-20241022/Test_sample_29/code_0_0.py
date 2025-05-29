# Create sets of friends
sarah_male_friends = {'Michael', 'Male2', 'Male3'}
sarah_female_friends = {'Female1', 'Female2', 'Female3'}

# For Michael, count his female friends
michael_female_friends = sarah_female_friends  # Michael is friends with all Sarah's female friends

print(f"Michael has {len(michael_female_friends)} female friends")