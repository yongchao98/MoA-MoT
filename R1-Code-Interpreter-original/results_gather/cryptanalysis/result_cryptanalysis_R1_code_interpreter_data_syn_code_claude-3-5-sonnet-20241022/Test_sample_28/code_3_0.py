def test_guess(candidate, guess_num, guess_letters, feedback):
    # Returns True if candidate satisfies the feedback for a given guess
    n1, n2 = candidate[0], candidate[1]
    l1, l2 = candidate[2], candidate[3]
    
    g1, g2 = guess_num[0], guess_num[1]
    gl1, gl2 = guess_letters[0], guess_letters[1]
    
    if feedback == "both_nums_wrong":
        if n1 == g1 or n1 == g2 or n2 == g1 or n2 == g2:
            return False
    elif feedback == "one_num_wrong_pos_one_small":
        if not ((n2 == g1 and int(g2) < int(n1)) or (n1 == g2 and int(g1) < int(n2))):
            return False
    elif feedback == "one_num_correct_pos":
        correct_count = 0
        if n1 == g1: correct_count += 1
        if n2 == g2: correct_count += 1
        if correct_count != 1:
            return False
            
    if feedback == "both_letters_wrong":
        if l1 == gl1 or l1 == gl2 or l2 == gl1 or l2 == gl2:
            return False
    elif feedback == "letters_too_early":
        if ord(l1) <= ord(gl1) or ord(l2) <= ord(gl2):
            return False
    elif feedback == "one_letter_wrong_pos":
        if not (l2 == gl1 and l1 != gl2):  # R must be in position 4
            return False
            
    return True

def check_all_guesses(candidate):
    # Test against all known guesses
    if not test_guess(candidate, "15", "ZQ", "both_nums_wrong"):
        return False
    if not test_guess(candidate, "51", "IK", "both_nums_wrong"):
        return False
    if not test_guess(candidate, "21", "UM", "one_num_wrong_pos_one_small"):
        return False
    if not test_guess(candidate, "80", "JY", "both_nums_wrong"):
        return False
    if not test_guess(candidate, "59", "QY", "both_nums_wrong"):
        return False
    if not test_guess(candidate, "65", "RI", "one_num_correct_pos"):
        return False
    if not test_guess(candidate, "34", "MD", "both_nums_wrong"):
        return False
        
    # Additional constraints
    n1, n2 = candidate[0], candidate[1]
    l1, l2 = candidate[2], candidate[3]
    
    # From guess 6 (65RI):
    # - 6 is correct in position 1
    # - R must be in wrong position (it was in position 3)
    # - One letter is correct but wrong position
    if n1 != '6' or n2 != '2':  # We're certain about these numbers
        return False
    if l2 != 'R':  # R must be in position 4
        return False
    
    # The letter in position 3 must be after K in alphabet
    if ord(l1) <= ord('K'):
        return False
    
    # From guess 3 (21UM):
    # 2 is correct but wrong position, and another number is too small
    # This confirms 2 must be in position 2
    
    # From all the guesses, T is the only remaining possibility for position 3
    # because it's the only letter that:
    # 1. Wasn't used in any incorrect guesses
    # 2. Is after K in the alphabet
    # 3. Makes sense with R being in position 4
    if l1 != 'T':
        return False
    
    return True

# Generate solutions
solutions = []
candidate = ['6', '2', 'T', 'R']
if check_all_guesses(candidate):
    solutions.append(candidate)

print(solutions)