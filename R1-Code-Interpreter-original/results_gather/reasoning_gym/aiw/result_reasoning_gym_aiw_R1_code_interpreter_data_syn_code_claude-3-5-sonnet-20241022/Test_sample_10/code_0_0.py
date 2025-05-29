# Let's create a simple representation of the friendship network
# M1 (William) and M2 are male friends
# F1 and F2 are female friends
# Mary is represented as 'Mary'

# Create a dictionary where each person has a list of their friends
friendships = {
    'Mary': ['M1', 'M2', 'F1', 'F2'],
    'M1': ['Mary', 'M2', 'F1', 'F2'],  # M1 is William
    'M2': ['Mary', 'M1', 'F1', 'F2'],
    'F1': ['Mary', 'M1', 'M2', 'F2'],
    'F2': ['Mary', 'M1', 'M2', 'F1']
}

# Count William's (M1's) female friends
william_female_friends = sum(1 for friend in friendships['M1'] if friend.startswith('F'))
print(f"William has {william_female_friends} female friends")