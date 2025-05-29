# Create sets of males and females
females = {'Linda', 'Female_Friend'}
males = {'Bob', 'Other_Male'}

# Create a function to count female friends of a person
def count_female_friends(person, females, males):
    # If the person is male, they are friends with all females
    # If the person is female, they are friends with all females except themselves
    if person in males:
        return len(females)
    else:
        return len(females) - 1

# Count Bob's female friends
bob_female_friends = count_female_friends('Bob', females, males)
print(f"Bob has {bob_female_friends} female friends")