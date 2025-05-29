from itertools import permutations

# Define the possible options
names = ['carol', 'alice', 'arnold', 'bob']
drinks = ['water', 'milk', 'boba tea', 'coffee']
phones = ['oneplus 9', 'samsung galaxy s21', 'huawei p50', 'google pixel 6']
pets = ['fish', 'bird', 'cat', 'dog']

# Define the constraints based on the clues
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (drink1, drink2, drink3, drink4) = arrangement['drinks']
    (phone1, phone2, phone3, phone4) = arrangement['phones']
    (pet1, pet2, pet3, pet4) = arrangement['pets']
    
    # Apply the constraints
    if name4 != 'bob' or drink4 != 'coffee':  # Clue 1
        return False
    if name3 != 'alice' or drink3 != 'milk':  # Clue 3 and 11
        return False
    if name1 != 'carol' or pet1 != 'bird':  # Clue 4 and 9
        return False
    if phone2 != 'samsung galaxy s21' or drink2 != 'water':  # Clue 8
        return False
    if phone3 != 'samsung galaxy s21' or drink4 != 'boba tea':  # Clue 7
        return False
    if phone4 != 'google pixel 6' or pet4 != 'fish':  # Clue 6
        return False
    if phone3 != 'oneplus 9' and phone4 != 'oneplus 9':  # Clue 2
        return False
    if pet3 != 'cat' and pet4 != 'cat':  # Clue 10
        return False
    if phone1 != 'google pixel 6' and phone2 != 'google pixel 6':  # Clue 5
        return False
    
    return True

# Try all permutations of the options
for name_perm in permutations(names):
    for drink_perm in permutations(drinks):
        for phone_perm in permutations(phones):
            for pet_perm in permutations(pets):
                arrangement = {
                    'names': name_perm,
                    'drinks': drink_perm,
                    'phones': phone_perm,
                    'pets': pet_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break