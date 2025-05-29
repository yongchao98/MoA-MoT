import string

def check_letter_order(l1, l2):
    # Check if letters satisfy all ordering conditions
    
    # From guess 1 (18VC): one correct position, one too early
    vc_match = 0
    vc_early = 0
    if l1 == 'V': vc_match += 1
    if l2 == 'C': vc_match += 1
    if l1 < 'V': vc_early += 1
    if l2 < 'C': vc_early += 1
    if vc_match != 1 or vc_early != 1:
        return False
    
    # From guess 2 (54OF): both too early
    if not (l1 > 'O' and l2 > 'F'):
        return False
        
    # From guess 3 (36AV): V in wrong position, A too early
    v_wrong_pos = False
    if l1 == 'V' or l2 == 'V':
        v_wrong_pos = True
    if not v_wrong_pos:
        return False
    if l1 < 'A' or l2 < 'A':
        return False
        
    # From guess 4 (42GE): both too early
    if not (l1 > 'G' and l2 > 'E'):
        return False
        
    # From guess 5 (08NQ): both incorrect
    if l1 in 'NQ' or l2 in 'NQ':
        return False
        
    return True

# Generate all possible combinations
valid_passwords = []
letters = string.ascii_uppercase

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # Letters can't repeat
            if check_letter_order(l1, l2):
                valid_passwords.append(['3', '6', l1, l2])

print(valid_passwords)