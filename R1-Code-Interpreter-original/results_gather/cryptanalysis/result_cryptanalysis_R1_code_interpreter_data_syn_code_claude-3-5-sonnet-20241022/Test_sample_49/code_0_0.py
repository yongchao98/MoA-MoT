def verify_guess(actual, guess, feedback):
    # Unpack the feedback
    num_correct = feedback['numbers_correct']
    num_wrong_pos = feedback['numbers_wrong_pos']
    num_too_large = feedback['numbers_too_large']
    num_too_small = feedback['numbers_too_small']
    letters_correct = feedback['letters_correct']
    
    # Check numbers
    correct_pos = 0
    wrong_pos = 0
    too_large = 0
    too_small = 0
    
    for i in range(2):
        if guess[i] == actual[i]:
            correct_pos += 1
        elif guess[i] in actual[:2]:
            wrong_pos += 1
        elif int(guess[i]) > int(actual[i]):
            too_large += 1
        else:
            too_small += 1
    
    # Check letters
    letters_match = sum(1 for i in range(2,4) if guess[i] == actual[i])
    
    # Verify all conditions
    conditions_met = (
        (num_correct == correct_pos) and
        (num_wrong_pos == wrong_pos) and
        (num_too_large == too_large) and
        (num_too_small == too_small) and
        (letters_correct == letters_match)
    )
    
    return conditions_met

# The proposed solution
solution = ['5', '4', 'U', 'A']

# All the guesses and their feedback
guesses = [
    {
        'guess': ['7', '5', 'K', 'B'],
        'feedback': {
            'numbers_correct': 0,
            'numbers_wrong_pos': 1,
            'numbers_too_large': 1,
            'numbers_too_small': 0,
            'letters_correct': 0
        }
    },
    {
        'guess': ['1', '5', 'U', 'A'],
        'feedback': {
            'numbers_correct': 0,
            'numbers_wrong_pos': 1,
            'numbers_too_small': 1,
            'numbers_too_large': 0,
            'letters_correct': 2
        }
    },
    {
        'guess': ['0', '9', 'K', 'J'],
        'feedback': {
            'numbers_correct': 0,
            'numbers_wrong_pos': 0,
            'numbers_too_large': 0,
            'numbers_too_small': 0,
            'letters_correct': 0
        }
    },
    {
        'guess': ['0', '5', 'P', 'B'],
        'feedback': {
            'numbers_correct': 0,
            'numbers_wrong_pos': 1,
            'numbers_too_small': 1,
            'numbers_too_large': 0,
            'letters_correct': 0
        }
    },
    {
        'guess': ['6', '9', 'D', 'B'],
        'feedback': {
            'numbers_correct': 0,
            'numbers_wrong_pos': 0,
            'numbers_too_large': 2,
            'numbers_too_small': 0,
            'letters_correct': 0
        }
    },
    {
        'guess': ['5', '4', 'Z', 'F'],
        'feedback': {
            'numbers_correct': 2,
            'numbers_wrong_pos': 0,
            'numbers_too_large': 0,
            'numbers_too_small': 0,
            'letters_correct': 0
        }
    }
]

# Verify the solution against all guesses
all_correct = all(verify_guess(solution, g['guess'], g['feedback']) for g in guesses)
print(f"Solution {solution} is {'correct' if all_correct else 'incorrect'} according to all conditions")