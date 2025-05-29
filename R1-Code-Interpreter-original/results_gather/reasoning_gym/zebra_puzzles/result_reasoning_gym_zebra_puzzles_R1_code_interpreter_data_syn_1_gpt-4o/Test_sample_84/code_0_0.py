from itertools import permutations

# Define the possible options for each category
names = ['arnold', 'carol', 'bob', 'alice']
phones = ['samsung galaxy s21', 'google pixel 6', 'oneplus 9', 'huawei p50']
pets = ['dog', 'cat', 'fish', 'bird']
drinks = ['coffee', 'boba tea', 'water', 'milk']

# Define the constraints based on the clues
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (phone1, phone2, phone3, phone4) = arrangement['phones']
    (pet1, pet2, pet3, pet4) = arrangement['pets']
    (drink1, drink2, drink3, drink4) = arrangement['drinks']
    
    # Apply the constraints
    if phone3 != 'oneplus 9' or drink3 != 'coffee':
        return False
    if drink1 != 'boba tea' or pet1 != 'fish':
        return False
    if name4 != 'alice' or pet4 != 'cat':
        return False
    if name4 != 'alice' or drink4 != 'water':
        return False
    if (drink2 != 'milk' or pet3 != 'bird') and (drink3 != 'milk' or pet2 != 'bird'):
        return False
    if name2 != 'carol' or phone2 != 'samsung galaxy s21':
        return False
    if name3 != 'arnold':
        return False
    if (phone2 != 'google pixel 6' or name1 != 'bob') and (phone4 != 'google pixel 6' or name3 != 'bob'):
        return False
    
    return True

# Generate all possible permutations of the arrangements
for names_perm in permutations(names):
    for phones_perm in permutations(phones):
        for pets_perm in permutations(pets):
            for drinks_perm in permutations(drinks):
                arrangement = {
                    'names': names_perm,
                    'phones': phones_perm,
                    'pets': pets_perm,
                    'drinks': drinks_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break