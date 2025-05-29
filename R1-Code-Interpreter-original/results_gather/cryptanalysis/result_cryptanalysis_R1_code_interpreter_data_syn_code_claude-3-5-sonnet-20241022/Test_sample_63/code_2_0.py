def check_numbers(n1, n2, guess_n1, guess_n2, correct_pos, wrong_pos, too_large, too_small):
    pos_correct = 0
    pos_wrong = 0
    larger = 0
    smaller = 0
    
    # Check first digit
    if n1 == guess_n1:
        pos_correct += 1
    elif n1 == guess_n2:
        pos_wrong += 1
    elif n1 > guess_n1:
        larger += 1
    else:
        smaller += 1
        
    # Check second digit
    if n2 == guess_n2:
        pos_correct += 1
    elif n2 == guess_n1:
        pos_wrong += 1
    elif n2 > guess_n2:
        larger += 1
    else:
        smaller += 1
        
    return (pos_correct == correct_pos and 
            pos_wrong == wrong_pos and 
            larger == too_large and 
            smaller == too_small)

def check_letters(l1, l2, guess_l1, guess_l2, correct_pos, wrong_pos, too_early, too_late):
    pos_correct = 0
    pos_wrong = 0
    early = 0
    late = 0
    
    # Check first letter
    if l1 == guess_l1:
        pos_correct += 1
    elif l1 == guess_l2:
        pos_wrong += 1
    elif ord(l1) < ord(guess_l1):
        early += 1
    else:
        late += 1
        
    # Check second letter
    if l2 == guess_l2:
        pos_correct += 1
    elif l2 == guess_l1:
        pos_wrong += 1
    elif ord(l2) < ord(guess_l2):
        early += 1
    else:
        late += 1
        
    return (pos_correct == correct_pos and 
            pos_wrong == wrong_pos and 
            early == too_early and 
            late == too_late)

# Store all guesses and their conditions
guesses = [
    ((3,8,'U','W'), (0,1,1,0), (0,0,0,2)),  # 38UW
    ((6,9,'O','G'), (1,0,1,0), (1,0,1,0)),  # 69OG
    ((8,3,'O','P'), (1,0,1,0), (1,0,0,0)),  # 83OP
    ((4,8,'S','P'), (0,0,0,0), (0,0,0,0)),  # 48SP
    ((4,2,'Z','E'), (0,0,0,0), (0,0,0,0)),  # 42ZE
    ((0,1,'W','G'), (0,0,0,2), (0,0,0,0)),  # 01WG
    ((0,7,'J','X'), (0,0,0,0), (0,0,0,0)),  # 07JX
    ((3,1,'H','K'), (0,1,0,1), (0,0,2,0)),  # 31HK
    ((3,7,'V','Y'), (0,1,1,0), (0,0,0,2)),  # 37VY
    ((6,2,'Q','H'), (1,0,0,1), (0,0,0,0)),  # 62QH
    ((9,4,'C','L'), (0,0,0,0), (0,0,2,0)),  # 94CL
    ((7,1,'T','E'), (0,0,0,0), (0,1,1,0)),  # 71TE
]

solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                if l1 == l2:
                    continue
                    
                valid = True
                for guess in guesses:
                    if not check_numbers(n1, n2, guess[0][0], guess[0][1], *guess[1]):
                        valid = False
                        break
                    if not check_letters(l1, l2, guess[0][2], guess[0][3], *guess[2]):
                        valid = False
                        break
                
                if valid:
                    solutions.append([n1, n2, l1, l2])

print("Found solutions:", solutions)

# Additional analysis of key patterns
print("\nAnalyzing key patterns:")
print("Guesses with 6 in correct position:", ["69OG", "62QH"])
print("Guesses with 3 in wrong position:", ["38UW", "31HK", "37VY"])
print("Guesses with O in correct position:", ["69OG", "83OP"])
print("Guesses with T involved:", ["71TE"])