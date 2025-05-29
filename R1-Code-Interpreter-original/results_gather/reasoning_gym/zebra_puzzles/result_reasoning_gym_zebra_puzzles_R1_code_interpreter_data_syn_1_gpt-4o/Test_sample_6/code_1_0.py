from itertools import permutations

# Define the possible options for each category
names = ['arnold', 'carol', 'bob', 'alice']
phones = ['huawei p50', 'samsung galaxy s21', 'oneplus 9', 'google pixel 6']
cigars = ['pall mall', 'dunhill', 'prince', 'blue master']
children = ['billy', 'timothy', 'bella', 'alice']

# Define the constraints based on the clues
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (phone1, phone2, phone3, phone4) = arrangement['phones']
    (cigar1, cigar2, cigar3, cigar4) = arrangement['cigars']
    (child1, child2, child3, child4) = arrangement['children']
    
    # Apply the constraints
    if name4 != 'alice' or phone4 != 'google pixel 6':
        return False
    if 'prince' not in [cigar1, cigar2, cigar3, cigar4] or 'huawei p50' not in [phone1, phone2, phone3, phone4]:
        return False
    if cigar3 != 'blue master' or name3 != 'bob' or phone3 != 'samsung galaxy s21':
        return False
    if name3 != 'bob' or child3 != 'billy':
        return False
    if name2 != 'carol' or child2 != 'alice':
        return False
    if cigar1 != 'pall mall' or child1 != 'alice':
        return False
    if name2 != 'carol' or cigar3 != 'prince':
        return False
    if name2 != 'carol' or name3 != 'bob':
        return False
    if name1 != 'arnold':
        return False
    
    return True

# Generate all possible permutations of the options
for names_perm in permutations(names):
    for phones_perm in permutations(phones):
        for cigars_perm in permutations(cigars):
            for children_perm in permutations(children):
                arrangement = {
                    'names': names_perm,
                    'phones': phones_perm,
                    'cigars': cigars_perm,
                    'children': children_perm
                }
                if is_valid(arrangement):
                    print(arrangement['names'][0])
                    break