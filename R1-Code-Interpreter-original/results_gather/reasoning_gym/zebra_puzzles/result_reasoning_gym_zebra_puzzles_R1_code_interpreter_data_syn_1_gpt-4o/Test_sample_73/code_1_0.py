from itertools import permutations

# Define the possible options
names = ['arnold', 'alice', 'bob', 'carol']
drinks = ['boba tea', 'milk', 'water', 'coffee']
flowers = ['daffodils', 'lilies', 'carnations', 'iris']
phones = ['oneplus 9', 'samsung galaxy s21', 'huawei p50', 'google pixel 6']

# Define the constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4), (drink1, drink2, drink3, drink4), (flower1, flower2, flower3, flower4), (phone1, phone2, phone3, phone4) = arrangement
    
    # Apply the constraints
    if name3 != 'carol' or drink3 != 'boba tea':
        return False
    if phone4 != 'oneplus 9':
        return False
    if flower1 != 'lilies':
        return False
    if name1 == 'arnold' and flower1 != 'carnations':
        return False
    if name2 == 'arnold' and flower2 != 'carnations':
        return False
    if name3 == 'arnold' and flower3 != 'carnations':
        return False
    if name4 == 'arnold' and flower4 != 'carnations':
        return False
    if name1 == 'bob' and drink1 != 'water':
        return False
    if name2 == 'bob' and drink2 != 'water':
        return False
    if name3 == 'bob' and drink3 != 'water':
        return False
    if name4 == 'bob' and drink4 != 'water':
        return False
    if phone1 == 'samsung galaxy s21' and flower2 != 'daffodils':
        return False
    if phone2 == 'samsung galaxy s21' and flower3 != 'daffodils':
        return False
    if phone3 == 'samsung galaxy s21' and flower4 != 'daffodils':
        return False
    if phone4 == 'samsung galaxy s21':
        return False
    if phone1 == 'google pixel 6' and drink1 != 'milk':
        return False
    if phone2 == 'google pixel 6' and drink2 != 'milk':
        return False
    if phone3 == 'google pixel 6' and drink3 != 'milk':
        return False
    if phone4 == 'google pixel 6' and drink4 != 'milk':
        return False
    if flower1 == 'lilies' and name2 != 'alice' and name2 != 'arnold':
        return False
    if flower2 == 'lilies' and name1 != 'alice' and name1 != 'arnold':
        return False
    if flower2 == 'lilies' and name3 != 'alice' and name3 != 'arnold':
        return False
    if flower3 == 'lilies' and name2 != 'alice' and name2 != 'arnold':
        return False
    if flower3 == 'lilies' and name4 != 'alice' and name4 != 'arnold':
        return False
    if flower4 == 'lilies' and name3 != 'alice' and name3 != 'arnold':
        return False
    return True

# Generate all possible permutations
for arrangement in permutations(permutations(names), 4):
    if is_valid(arrangement):
        print(arrangement)
        break