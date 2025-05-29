def is_valid_numbers(n1, n2, guess, feedback):
    if feedback == "both_small":
        return n1 > int(guess[0]) and n2 > int(guess[1])
    elif feedback == "one_correct_wrong_pos_one_small":
        # One number is correct but in wrong position
        correct_wrong_pos = (str(n1) == guess[1] or str(n2) == guess[0])
        # One number is too small
        one_small = (n1 < int(guess[0]) or n2 < int(guess[1]))
        return correct_wrong_pos and one_small
    elif feedback == "both_incorrect":
        # Neither number should appear in the guess
        return str(n1) not in guess and str(n2) not in guess
    return False

def is_valid_letters(l1, l2, guess, feedback):
    if feedback == "one_correct_wrong_pos_one_early":
        # One letter is correct but in wrong position
        correct_wrong_pos = (l1 == guess[1] or l2 == guess[0])
        # One letter is too early
        one_early = (ord(l1) < ord(guess[0]) or ord(l2) < ord(guess[1]))
        return correct_wrong_pos and one_early
    elif feedback == "both_early":
        return ord(l1) > ord(guess[0]) and ord(l2) > ord(guess[1])
    elif feedback == "both_incorrect":
        return l1 not in guess and l2 not in guess
    return False

def find_password():
    # Try all possible number combinations
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:  # Numbers can't repeat
                continue
                
            # Check all number conditions
            if not all([
                is_valid_numbers(n1, n2, "30", "both_small"),
                is_valid_numbers(n1, n2, "60", "one_correct_wrong_pos_one_small"),
                is_valid_numbers(n1, n2, "23", "both_small"),
                is_valid_numbers(n1, n2, "03", "both_small"),
                is_valid_numbers(n1, n2, "79", "both_incorrect")
            ]):
                continue

            # Try all possible letter combinations
            for l1 in [chr(i) for i in range(ord('A'), ord('Z')+1)]:
                for l2 in [chr(i) for i in range(ord('A'), ord('Z')+1)]:
                    if l1 == l2:  # Letters can't repeat
                        continue
                        
                    # Check all letter conditions
                    if all([
                        is_valid_letters(l1, l2, "AU", "one_correct_wrong_pos_one_early"),
                        is_valid_letters(l1, l2, "KI", "both_early"),
                        is_valid_letters(l1, l2, "QB", "both_incorrect"),
                        is_valid_letters(l1, l2, "HC", "both_early"),
                        is_valid_letters(l1, l2, "MI", "both_incorrect")
                    ]):
                        return [str(n1), str(n2), l1, l2]
    return None

result = find_password()
print(result)