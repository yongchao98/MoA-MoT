# Define the people
people = {
    'Alice': 'F',
    'Richard': 'M',
    'Male2': 'M',  # The other male friend
    'Female1': 'F'  # Alice's female friend
}

# Count Richard's female friends
richard_female_friends = sum(1 for person, gender in people.items() 
                           if gender == 'F' and person != 'Richard')

print(f"Richard has {richard_female_friends} female friends")