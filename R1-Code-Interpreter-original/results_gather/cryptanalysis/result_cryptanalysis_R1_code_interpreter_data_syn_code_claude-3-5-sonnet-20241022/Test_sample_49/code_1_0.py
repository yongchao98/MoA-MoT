def verify_guess(actual, guess, feedback):
    # Convert numbers to integers for comparison
    actual_nums = [int(actual[0]), int(actual[1])]
    guess_nums = [int(guess[0]), int(guess[1])]
    
    # Count correct numbers in correct positions
    correct_nums = sum(1 for i in range(2) if guess_nums[i] == actual_nums[i])
    
    # Count numbers in wrong positions
    wrong_pos = sum(1 for i in range(2) for j in range(2) if i != j and guess_nums[i] == actual_nums[j])
    
    # Count too large/small numbers
    too_large = sum(1 for i in range(2) if guess_nums[i] > actual_nums[i] and guess_nums[i] not in actual_nums)
    too_small = sum(1 for i in range(2) if guess_nums[i] < actual_nums[i] and guess_nums[i] not in actual_nums)
    
    # Count correct letters
    correct_letters = sum(1 for i in range(2,4) if guess[i] == actual[i])
    
    # Check if feedback matches
    return (
        correct_nums == feedback['correct_nums'] and
        wrong_pos == feedback['wrong_pos'] and
        too_large == feedback['too_large'] and
        too_small == feedback['too_small'] and
        correct_letters == feedback['correct_letters']
    )

def try_all_combinations():
    numbers = ['5', '4']  # We know these are correct from guess 6
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    feedbacks = [
        # 75KB
        {'guess': ['7', '5', 'K', 'B'], 
         'feedback': {'correct_nums': 0, 'wrong_pos': 1, 'too_large': 1, 'too_small': 0, 'correct_letters': 0}},
        # 15UA
        {'guess': ['1', '5', 'U', 'A'], 
         'feedback': {'correct_nums': 0, 'wrong_pos': 1, 'too_large': 0, 'too_small': 1, 'correct_letters': 2}},
        # 09KJ
        {'guess': ['0', '9', 'K', 'J'], 
         'feedback': {'correct_nums': 0, 'wrong_pos': 0, 'too_large': 0, 'too_small': 2, 'correct_letters': 0}},
        # 05PB
        {'guess': ['0', '5', 'P', 'B'], 
         'feedback': {'correct_nums': 0, 'wrong_pos': 1, 'too_large': 0, 'too_small': 1, 'correct_letters': 0}},
        # 69DB
        {'guess': ['6', '9', 'D', 'B'], 
         'feedback': {'correct_nums': 0, 'wrong_pos': 0, 'too_large': 2, 'too_small': 0, 'correct_letters': 0}},
        # 54ZF
        {'guess': ['5', '4', 'Z', 'F'], 
         'feedback': {'correct_nums': 2, 'wrong_pos': 0, 'too_large': 0, 'too_small': 0, 'correct_letters': 0}}
    ]
    
    for l1 in letters:
        for l2 in letters:
            if l1 != l2:  # Letters can't repeat
                combination = ['5', '4', l1, l2]
                if all(verify_guess(combination, f['guess'], f['feedback']) for f in feedbacks):
                    print(f"Found valid combination: {combination}")

try_all_combinations()