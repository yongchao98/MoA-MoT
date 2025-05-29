from itertools import permutations

# Define the possible options
names = ['alice', 'arnold', 'bob', 'carol']
cigars = ['prince', 'dunhill', 'blue master', 'pall mall']
drinks = ['milk', 'boba tea', 'coffee', 'water']
phones = ['huawei p50', 'oneplus 9', 'google pixel 6', 'samsung galaxy s21']

# Define the constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (cigar1, cigar2, cigar3, cigar4) = arrangement['cigars']
    (drink1, drink2, drink3, drink4) = arrangement['drinks']
    (phone1, phone2, phone3, phone4) = arrangement['phones']
    
    # Apply the constraints
    if cigar1 != 'dunhill': return False
    if phone4 != 'samsung galaxy s21' or drink4 != 'coffee': return False
    if name3 != 'alice' or cigar3 != 'blue master' or drink3 != 'milk': return False
    if name4 != 'carol' or drink4 != 'water': return False
    if name2 != 'bob' or phone2 != 'oneplus 9': return False
    if (cigar3 == 'blue master' and phone2 != 'oneplus 9') and (cigar2 == 'blue master' and phone3 != 'oneplus 9'): return False
    if (drink2 == 'boba tea' and cigar3 != 'prince') and (drink3 == 'boba tea' and cigar2 != 'prince'): return False
    if (drink1 == 'milk' and phone1 != 'huawei p50') and (drink2 == 'milk' and phone2 != 'huawei p50') and (drink3 == 'milk' and phone3 != 'huawei p50') and (drink4 == 'milk' and phone4 != 'huawei p50'): return False
    
    return True

# Generate all permutations and check for valid arrangements
for name_perm in permutations(names):
    for cigar_perm in permutations(cigars):
        for drink_perm in permutations(drinks):
            for phone_perm in permutations(phones):
                arrangement = {
                    'names': name_perm,
                    'cigars': cigar_perm,
                    'drinks': drink_perm,
                    'phones': phone_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break