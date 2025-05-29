def check_guess(password, guess, feedback):
    # Unpack the guess and feedback
    n1, n2 = int(guess[0]), int(guess[1])
    l1, l2 = guess[2], guess[3]
    pn1, pn2 = int(password[0]), int(password[1])
    pl1, pl2 = password[2], password[3]
    
    if feedback == "both_numbers_wrong":
        if n1 in [pn1, pn2] or n2 in [pn1, pn2]:
            return False
    elif feedback == "one_number_right_wrong_pos_one_too_large":
        count_correct = (n1 == pn2) + (n2 == pn1)
        if count_correct != 1:
            return False
        if n1 > pn1 and n1 > pn2 and n2 > pn1 and n2 > pn2:
            return False
    elif feedback == "both_numbers_right_wrong_pos":
        if not (n1 == pn2 and n2 == pn1):
            return False
    
    if feedback == "one_letter_right_wrong_pos_one_wrong":
        count_correct = (l1 == pl2) + (l2 == pl1)
        if count_correct != 1:
            return False
    elif feedback == "both_letters_wrong":
        if l1 in [pl1, pl2] or l2 in [pl1, pl2]:
            return False
    elif feedback == "both_letters_too_late":
        if l1 <= pl1 and l1 <= pl2 and l2 <= pl1 and l2 <= pl2:
            return False
    elif feedback == "one_letter_right_pos_one_too_late":
        if not ((l1 == pl1 and l2 > max(pl1, pl2)) or (l2 == pl2 and l1 > max(pl1, pl2))):
            return False
    
    return True

def find_password():
    numbers = range(0, 10)
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    password = [str(n1), str(n2), l1, l2]
                    
                    # Check all clues
                    if not check_guess(password, "14HF", "both_numbers_wrong"):
                        continue
                    if not check_guess(password, "38FX", "one_number_right_wrong_pos_one_too_large"):
                        continue
                    if not check_guess(password, "52YF", "one_number_right_wrong_pos_one_too_large"):
                        continue
                    if not check_guess(password, "32QI", "both_numbers_right_wrong_pos"):
                        continue
                    if not check_guess(password, "71DX", "both_numbers_wrong"):
                        continue
                    
                    # Additional letter checks
                    if 'H' in [l1, l2] and 'F' not in [l1, l2]:  # From clue 1
                        continue
                    if 'F' in [l1, l2]:  # From clue 2
                        continue
                    if 'Y' in [l1, l2]:  # From clue 3
                        continue
                    if 'Q' in [l1, l2] or 'I' in [l1, l2]:  # From clue 4
                        continue
                    if 'X' in [l1, l2]:  # From clue 5
                        continue
                    
                    print(password)

find_password()