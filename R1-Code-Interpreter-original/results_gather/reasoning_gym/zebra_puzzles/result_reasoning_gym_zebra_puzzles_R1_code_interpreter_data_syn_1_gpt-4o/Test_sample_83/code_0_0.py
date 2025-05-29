from itertools import permutations

names = ['arnold', 'carol', 'bob', 'alice']
lunches = ['stir fry', 'pizza', 'soup', 'grilled cheese']
phones = ['oneplus 9', 'samsung galaxy s21', 'google pixel 6', 'huawei p50']
children = ['timothy', 'bella', 'alice', 'billy']

# Constraints
def is_valid(arrangement):
    name, lunch, phone, child = arrangement
    
    # Bob uses Huawei P50
    if name.index('bob') != phone.index('huawei p50'):
        return False
    
    # Alice loves stir fry and is the mother of Billy
    if name.index('alice') != lunch.index('stir fry') or name.index('alice') != child.index('billy'):
        return False
    
    # Bob loves grilled cheese
    if name.index('bob') != lunch.index('grilled cheese'):
        return False
    
    # The person who loves pizza uses OnePlus 9 and is the mother of Timothy
    if lunch.index('pizza') != phone.index('oneplus 9') or lunch.index('pizza') != child.index('timothy'):
        return False
    
    # Carol is directly left of the person who uses a Samsung Galaxy S21
    if not (0 <= name.index('carol') < 3 and name.index('carol') + 1 == phone.index('samsung galaxy s21')):
        return False
    
    # Bob is directly left of the mother of Timothy
    if not (0 <= name.index('bob') < 3 and name.index('bob') + 1 == child.index('timothy')):
        return False
    
    # The person who loves soup and the person whose child is Bella are next to each other
    if abs(lunch.index('soup') - child.index('bella')) != 1:
        return False
    
    return True

# Find the valid arrangement
for name_perm in permutations(names):
    for lunch_perm in permutations(lunches):
        for phone_perm in permutations(phones):
            for child_perm in permutations(children):
                arrangement = (name_perm, lunch_perm, phone_perm, child_perm)
                if is_valid(arrangement):
                    print("Valid arrangement found:")
                    print("Names:", name_perm)
                    print("Lunches:", lunch_perm)
                    print("Phones:", phone_perm)
                    print("Children:", child_perm)
                    break