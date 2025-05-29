# Create sets of males and females
females = {'Susan', 'Female_Friend1'}
males = {'Michael', 'Male_Friend2', 'Male_Friend3'}

# Count Michael's female friends
michaels_female_friends = set()

# Add all females except Michael himself
for person in females:
    if person != 'Michael':  # Not needed as Michael is not in females, but good practice
        michaels_female_friends.add(person)

print(f"Michael has {len(michaels_female_friends)} female friends")