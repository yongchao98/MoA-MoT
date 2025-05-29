# Let's create a dictionary to represent friendships
# Let's name the people:
# Linda (L)
# Charles (C) - male friend
# M2 - other male friend
# F1, F2 - female friends

friendships = {
    'Linda': ['Charles', 'M2', 'F1', 'F2'],
    'Charles': ['Linda', 'M2', 'F1', 'F2'],  # Charles is friends with everyone
    'M2': ['Linda', 'Charles', 'F1', 'F2'],
    'F1': ['Linda', 'Charles', 'M2', 'F2'],
    'F2': ['Linda', 'Charles', 'M2', 'F1']
}

# Count Charles's female friends
charles_female_friends = sum(1 for friend in friendships['Charles'] 
                           if friend.startswith('F'))

print(f"Charles has {charles_female_friends} female friends")