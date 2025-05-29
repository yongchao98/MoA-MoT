from itertools import product

def get_possible_numbers():
    all_numbers = set(range(10))
    
    # From conditions where both numbers are incorrect
    impossible_pairs = {(8,0), (2,5), (4,8), (8,0), (4,3), (0,3), (2,3), (0,5), (6,0)}
    
    # From condition 2: 91PE - one correct in position, one too large
    # From condition 3: 51DH - one correct in position
    # From condition 4: 19HC - one correct but wrong position
    # From condition 7: 21GK - one correct in position
    
    possible_first = set()
    possible_second = set()
    
    # Test all possible number combinations
    for n1, n2 in product(range(10), range(10)):
        if n1 == n2:  # Numbers can't repeat
            continue
        if (n1, n2) in impossible_pairs:
            continue
            
        # Apply conditions
        if n1 == 9 or n1 == 8 or n1 > 7:  # Various conditions indicate these are too large
            continue
        if n2 == 0 or n2 == 8:  # These are definitely wrong from multiple conditions
            continue
            
        possible_first.add(n1)
        possible_second.add(n2)
    
    return possible_first, possible_second

def get_possible_letters():
    all_letters = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    
    # From conditions where both letters are incorrect
    impossible_pairs = {('Z','V'), ('P','E'), ('D','H'), ('F','G'), ('G','K'),
                       ('F','L'), ('D','L'), ('Z','O'), ('T','Q'), ('Q','F')}
    
    # From condition 4: HC - C is correct in position
    # From condition 11: SC - C is correct in position
    # From condition 12: XO - O is correct in position
    
    possible_first = set()
    possible_second = set()
    
    # Test all possible letter combinations
    for l1, l2 in product(all_letters, all_letters):
        if l1 == l2:  # Letters can't repeat
            continue
        if (l1, l2) in impossible_pairs:
            continue
            
        # Apply specific conditions
        if l2 != 'C':  # Second letter must be C from multiple conditions
            continue
        
        possible_first.add(l1)
        possible_second.add(l2)
    
    return possible_first, possible_second

# Get possible values
num1, num2 = get_possible_numbers()
let1, let2 = get_possible_letters()

print("Possible first numbers:", sorted(list(num1)))
print("Possible second numbers:", sorted(list(num2)))
print("Possible first letters:", sorted(list(let1)))
print("Possible second letters:", sorted(list(let2)))