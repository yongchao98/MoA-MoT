# Create sets of friends
male_friends = 2
female_friends = 4

# Create a dictionary to represent John's connections
johns_friends = {
    'female_friends': female_friends,  # All of Sarah's female friends
    'male_friends': male_friends - 1   # All male friends except himself
}

print(f"John has {johns_friends['female_friends']} female friends")