# Create sets of friends
male_friends = {'James', 'Male2'}  # 2 male friends
female_friends = {'Female1', 'Female2', 'Female3'}  # 3 female friends

# For James (who we know is one of the male friends)
james_female_friends = female_friends  # James is friends with all females

print(f"James has {len(james_female_friends)} female friends")