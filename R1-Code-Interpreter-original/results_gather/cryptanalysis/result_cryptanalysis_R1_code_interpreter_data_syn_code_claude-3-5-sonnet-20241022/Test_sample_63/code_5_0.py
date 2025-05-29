def verify_guess(test, guess, num_feedback, letter_feedback):
    test_nums = [test[0], test[1]]
    test_letters = [test[2], test[3]]
    guess_nums = [int(guess[0]), int(guess[1])]
    guess_letters = [guess[2], guess[3]]
    
    # Number verification
    num_correct_pos = sum(1 for i in range(2) if test_nums[i] == guess_nums[i])
    num_wrong_pos = sum(1 for i in range(2) if test_nums[i] in guess_nums and test_nums[i] != guess_nums[i])
    num_too_large = sum(1 for i in range(2) if test_nums[i] > guess_nums[i] and test_nums[i] != guess_nums[i])
    num_too_small = sum(1 for i in range(2) if test_nums[i] < guess_nums[i] and test_nums[i] != guess_nums[i])
    
    # Letter verification
    letter_correct_pos = sum(1 for i in range(2) if test_letters[i] == guess_letters[i])
    letter_wrong_pos = sum(1 for i in range(2) if test_letters[i] in guess_letters and test_letters[i] != guess_letters[i])
    letter_too_early = sum(1 for i in range(2) if test_letters[i] < guess_letters[i] and test_letters[i] != guess_letters[i])
    letter_too_late = sum(1 for i in range(2) if test_letters[i] > guess_letters[i] and test_letters[i] != guess_letters[i])
    
    return (num_correct_pos == num_feedback[0] and 
            num_wrong_pos == num_feedback[1] and 
            num_too_large == num_feedback[2] and 
            num_too_small == num_feedback[3] and 
            letter_correct_pos == letter_feedback[0] and 
            letter_wrong_pos == letter_feedback[1] and 
            letter_too_early == letter_feedback[2] and 
            letter_too_late == letter_feedback[3])

# Test cases with their feedback
test_cases = [
    ('38UW', (0,1,1,0), (0,0,0,2)),
    ('69OG', (1,0,1,0), (1,0,1,0)),
    ('83OP', (1,0,1,0), (1,0,0,0)),
    ('48SP', (0,0,0,0), (0,0,0,0)),
    ('42ZE', (0,0,0,0), (0,0,0,0)),
    ('01WG', (0,0,0,2), (0,0,0,0)),
    ('07JX', (0,0,0,0), (0,0,0,0)),
    ('31HK', (0,1,0,1), (0,0,2,0)),
    ('37VY', (0,1,1,0), (0,0,0,2)),
    ('62QH', (1,0,0,1), (0,0,0,0)),
    ('94CL', (0,0,0,0), (0,0,2,0)),
    ('71TE', (0,0,0,0), (0,1,1,0))
]

# Test specific combinations based on our analysis
test_combinations = [
    [6,3,'T','O'],  # Our primary guess
    [6,3,'O','T'],  # Alternative arrangement
    [3,6,'T','O'],  # Alternative arrangement
    [3,6,'O','T']   # Alternative arrangement
]

for combo in test_combinations:
    print(f"\nTesting combination: {combo}")
    all_valid = True
    for case in test_cases:
        is_valid = verify_guess(combo, case[0], case[1], case[2])
        print(f"Case {case[0]}: {'Valid' if is_valid else 'Invalid'}")
        if not is_valid:
            all_valid = False
    if all_valid:
        print(f"\nFound valid solution: {combo}")

# Additional systematic search
for n1 in [3,6]:  # Most likely first digits based on clues
    for n2 in [3,6]:  # Most likely second digits
        if n1 == n2:
            continue
        for l1 in ['O','T']:  # Most likely letters based on clues
            for l2 in ['O','T']:
                if l1 == l2:
                    continue
                combo = [n1,n2,l1,l2]
                all_valid = True
                for case in test_cases:
                    if not verify_guess(combo, case[0], case[1], case[2]):
                        all_valid = False
                        break
                if all_valid:
                    print(f"\nFound valid solution through systematic search: {combo}")