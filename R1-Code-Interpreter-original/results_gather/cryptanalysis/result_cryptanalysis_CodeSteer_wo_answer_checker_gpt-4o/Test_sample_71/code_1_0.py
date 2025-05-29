from itertools import permutations

# Define possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Define the feedback conditions
conditions = [
    ('65HO', (1, 1, 0, 0)),  # one number correct but wrong position, one number too small, both letters incorrect
    ('47IM', (0, 0, 0, 0)),  # both numbers and letters incorrect
    ('10GI', (0, 0, 0, 0)),  # both numbers too small, both letters too early
    ('70XG', (0, 0, 0, 0)),  # both numbers and letters incorrect
    ('12JU', (0, 0, 1, 0)),  # both numbers too small, one letter correct and correct position, one letter too late
    ('52SZ', (0, 0, 0, 0)),  # both numbers too small, both letters too late
    ('64AV', (1, 1, 0, 0)),  # one number correct but wrong position, one number too small, both letters incorrect
    ('85BO', (1, 0, 0, 0)),  # one number correct and correct position, one number too small, both letters incorrect
    ('98RX', (1, 1, 0, 0)),  # one number correct but wrong position, one number too large, both letters too late
    ('60JU', (1, 1, 1, 0)),  # one number correct but wrong position, one number too small, one letter correct and correct position, one letter too late
    ('75UC', (0, 0, 0, 0)),  # both numbers and letters incorrect
    ('79NC', (0, 0, 1, 0)),  # both numbers incorrect, one letter correct but wrong position, one letter too early
]

# Function to check if a combination satisfies a condition
def satisfies_condition(combo, guess, feedback):
    num_correct_pos, num_wrong_pos, letter_correct_pos, letter_wrong_pos = feedback
    num1, num2, letter1, letter2 = combo
    g_num1, g_num2, g_letter1, g_letter2 = guess

    # Check numbers
    num_correct = (num1 == g_num1) + (num2 == g_num2)
    num_wrong = (num1 == g_num2) + (num2 == g_num1)
    if num_correct != num_correct_pos or num_wrong != num_wrong_pos:
        return False

    # Check letters
    letter_correct = (letter1 == g_letter1) + (letter2 == g_letter2)
    letter_wrong = (letter1 == g_letter2) + (letter2 == g_letter1)
    if letter_correct != letter_correct_pos or letter_wrong != letter_wrong_pos:
        return False

    return True

# Iterate over all possible combinations
for num_combo in permutations(numbers, 2):
    for letter_combo in permutations(letters, 2):
        combo = num_combo + letter_combo
        if all(satisfies_condition(combo, guess, feedback) for guess, feedback in conditions):
            print(f"<<< {list(combo)} >>>")
            break