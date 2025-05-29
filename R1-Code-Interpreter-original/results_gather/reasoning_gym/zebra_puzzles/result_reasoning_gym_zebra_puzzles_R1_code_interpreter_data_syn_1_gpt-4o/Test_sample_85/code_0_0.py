from itertools import permutations

# Define the possible options for each category
names = ['carol', 'bob', 'arnold', 'alice']
cigars = ['blue master', 'prince', 'pall mall', 'dunhill']
colors = ['brown', 'purple', 'white', 'blue']
phones = ['samsung galaxy s21', 'google pixel 6', 'huawei p50', 'oneplus 9']

# Define the constraints based on the clues
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (cigar1, cigar2, cigar3, cigar4) = arrangement['cigars']
    (color1, color2, color3, color4) = arrangement['colors']
    (phone1, phone2, phone3, phone4) = arrangement['phones']
    
    # Apply the constraints
    if name4 != 'arnold': return False
    if cigar4 != 'prince': return False
    if color4 != 'blue': return False
    if color1 != 'white': return False
    if color2 != 'purple': return False
    if name1 == 'alice' and cigar1 != 'blue master': return False
    if name2 == 'alice' and cigar2 != 'blue master': return False
    if name3 == 'alice' and cigar3 != 'blue master': return False
    if name4 == 'alice' and cigar4 != 'blue master': return False
    if phone2 != 'huawei p50': return False
    if (cigar1 == 'pall mall' and phone2 != 'google pixel 6') and (cigar2 == 'pall mall' and phone1 != 'google pixel 6') and (cigar2 == 'pall mall' and phone3 != 'google pixel 6') and (cigar3 == 'pall mall' and phone2 != 'google pixel 6') and (cigar3 == 'pall mall' and phone4 != 'google pixel 6') and (cigar4 == 'pall mall' and phone3 != 'google pixel 6'): return False
    if (phone1 == 'samsung galaxy s21' and name2 != 'bob') and (phone2 == 'samsung galaxy s21' and name1 != 'bob' and name3 != 'bob') and (phone3 == 'samsung galaxy s21' and name2 != 'bob' and name4 != 'bob') and (phone4 == 'samsung galaxy s21' and name3 != 'bob'): return False
    if (cigar1 == 'pall mall' and phone1 != 'oneplus 9') and (cigar2 == 'pall mall' and phone2 != 'oneplus 9') and (cigar3 == 'pall mall' and phone3 != 'oneplus 9') and (cigar4 == 'pall mall' and phone4 != 'oneplus 9'): return False
    
    return True

# Try all permutations of the arrangements
for name_perm in permutations(names):
    for cigar_perm in permutations(cigars):
        for color_perm in permutations(colors):
            for phone_perm in permutations(phones):
                arrangement = {
                    'names': name_perm,
                    'cigars': cigar_perm,
                    'colors': color_perm,
                    'phones': phone_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break