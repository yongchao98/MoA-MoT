from itertools import permutations

names = ['carol', 'arnold', 'alice', 'bob']
phones = ['huawei p50', 'samsung galaxy s21', 'oneplus 9', 'google pixel 6']
drinks = ['milk', 'boba tea', 'coffee', 'water']
animals = ['bird', 'cat', 'fish', 'dog']

# Constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (phone1, phone2, phone3, phone4) = arrangement['phones']
    (drink1, drink2, drink3, drink4) = arrangement['drinks']
    (animal1, animal2, animal3, animal4) = arrangement['animals']
    
    # Apply the constraints
    if name2 != 'alice' or animal2 != 'cat':
        return False
    if drink3 != 'milk':
        return False
    if name1 == 'bob' and phone1 != 'huawei p50':
        return False
    if name2 == 'bob' and phone2 != 'huawei p50':
        return False
    if name3 == 'bob' and phone3 != 'huawei p50':
        return False
    if name4 == 'bob' and phone4 != 'huawei p50':
        return False
    if drink1 == 'water' and animal1 != 'bird':
        return False
    if drink2 == 'water' and animal2 != 'bird':
        return False
    if drink3 == 'water' and animal3 != 'bird':
        return False
    if drink4 == 'water' and animal4 != 'bird':
        return False
    if drink1 == 'boba tea' and animal1 != 'dog':
        return False
    if drink2 == 'boba tea' and animal2 != 'dog':
        return False
    if drink3 == 'boba tea' and animal3 != 'dog':
        return False
    if drink4 == 'boba tea' and animal4 != 'dog':
        return False
    if phone1 == 'google pixel 6' and name2 != 'carol':
        return False
    if phone2 == 'google pixel 6' and name3 != 'carol':
        return False
    if phone3 == 'google pixel 6' and name4 != 'carol':
        return False
    if name1 == 'carol' and drink1 != 'water':
        return False
    if name2 == 'carol' and drink2 != 'water':
        return False
    if name3 == 'carol' and drink3 != 'water':
        return False
    if name4 == 'carol' and drink4 != 'water':
        return False
    if name1 == 'carol' and phone1 != 'oneplus 9':
        return False
    if name2 == 'carol' and phone2 != 'oneplus 9':
        return False
    if name3 == 'carol' and phone3 != 'oneplus 9':
        return False
    if name4 == 'carol' and phone4 != 'oneplus 9':
        return False
    
    return True

# Generate all permutations and check constraints
for name_perm in permutations(names):
    for phone_perm in permutations(phones):
        for drink_perm in permutations(drinks):
            for animal_perm in permutations(animals):
                arrangement = {
                    'names': name_perm,
                    'phones': phone_perm,
                    'drinks': drink_perm,
                    'animals': animal_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break