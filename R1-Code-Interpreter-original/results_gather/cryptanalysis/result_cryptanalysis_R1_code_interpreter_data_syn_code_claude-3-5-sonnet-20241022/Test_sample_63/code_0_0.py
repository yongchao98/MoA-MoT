def check_guess(guess, correct_pos_nums, wrong_pos_nums, too_large_nums, too_small_nums, 
                correct_pos_letters, wrong_pos_letters, too_early_letters, too_late_letters):
    return {
        'guess': guess,
        'correct_pos_nums': correct_pos_nums,
        'wrong_pos_nums': wrong_pos_nums,
        'too_large_nums': too_large_nums,
        'too_small_nums': too_small_nums,
        'correct_pos_letters': correct_pos_letters,
        'wrong_pos_letters': wrong_pos_letters,
        'too_early_letters': too_early_letters,
        'too_late_letters': too_late_letters
    }

guesses = [
    check_guess('38UW', 0, 1, 1, 0, 0, 0, 0, 2),
    check_guess('69OG', 1, 0, 1, 0, 1, 0, 1, 0),
    check_guess('83OP', 1, 0, 1, 0, 1, 0, 0, 0),
    check_guess('48SP', 0, 0, 0, 0, 0, 0, 0, 0),
    check_guess('42ZE', 0, 0, 0, 0, 0, 0, 0, 0),
    check_guess('01WG', 0, 0, 0, 2, 0, 0, 0, 0),
    check_guess('07JX', 0, 0, 0, 0, 0, 0, 0, 0),
    check_guess('31HK', 0, 1, 0, 1, 0, 0, 2, 0),
    check_guess('37VY', 0, 1, 1, 0, 0, 0, 0, 2),
    check_guess('62QH', 1, 0, 0, 1, 0, 0, 0, 0),
    check_guess('94CL', 0, 0, 0, 0, 0, 0, 2, 0),
    check_guess('71TE', 0, 0, 0, 0, 0, 1, 1, 0)
]

# Let's analyze the numbers first
numbers = set(range(10))
possible_numbers = []

for i in range(10):
    for j in range(10):
        if i != j:  # numbers can't repeat
            valid = True
            for g in guesses:
                nums = [int(g['guess'][0]), int(g['guess'][1])]
                correct_count = 0
                correct_pos = 0
                too_large = 0
                too_small = 0
                
                for pos in range(2):
                    if [i,j][pos] == nums[pos]:
                        correct_pos += 1
                    elif [i,j][pos] in nums:
                        correct_count += 1
                    elif [i,j][pos] > nums[pos]:
                        too_large += 1
                    else:
                        too_small += 1
                
                if (correct_pos != g['correct_pos_nums'] or 
                    correct_count != g['wrong_pos_nums'] or
                    too_large != g['too_large_nums'] or
                    too_small != g['too_small_nums']):
                    valid = False
                    break
            
            if valid:
                possible_numbers.append([i,j])

print("Possible number combinations:", possible_numbers)

# Now analyze letters
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
possible_letters = []

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # letters can't repeat
            valid = True
            for g in guesses:
                guess_letters = [g['guess'][2], g['guess'][3]]
                correct_pos = 0
                correct_wrong_pos = 0
                too_early = 0
                too_late = 0
                
                for pos in range(2):
                    if [l1,l2][pos] == guess_letters[pos]:
                        correct_pos += 1
                    elif [l1,l2][pos] in guess_letters:
                        correct_wrong_pos += 1
                    elif ord([l1,l2][pos]) < ord(guess_letters[pos]):
                        too_early += 1
                    else:
                        too_late += 1
                
                if (correct_pos != g['correct_pos_letters'] or 
                    correct_wrong_pos != g['wrong_pos_letters'] or
                    too_early != g['too_early_letters'] or
                    too_late != g['too_late_letters']):
                    valid = False
                    break
            
            if valid:
                possible_letters.append([l1,l2])

print("Possible letter combinations:", possible_letters)