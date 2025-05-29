# Create sets of friends
male_friends = set(['David', 'M2', 'M3', 'M4', 'M5'])  # Including David and other males
female_friends = set(['F1', 'F2', 'F3', 'F4', 'F5'])  # Female friends
jessica = 'Jessica'

# Count David's female friends
david_female_friends = female_friends  # Since everyone is friends with everyone

print(f"Number of David's female friends: {len(david_female_friends)}")