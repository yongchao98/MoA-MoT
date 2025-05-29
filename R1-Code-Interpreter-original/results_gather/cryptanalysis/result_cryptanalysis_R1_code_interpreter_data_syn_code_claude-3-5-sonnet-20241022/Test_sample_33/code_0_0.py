from itertools import permutations

def check_guess(guess, actual):
    # Check numbers
    num_correct_pos = 0
    num_wrong_pos = 0
    num_too_small = 0
    for i in range(2):
        if guess[i] == actual[i]:
            num_correct_pos += 1
        elif guess[i] in actual[:2]:
            num_wrong_pos += 1
        elif int(guess[i]) < int(actual[i]):
            num_too_small += 1
    
    # Check letters
    let_correct_pos = 0
    let_wrong_pos = 0
    for i in range(2,4):
        if guess[i] == actual[i]:
            let_correct_pos += 1
        elif guess[i] in actual[2:]:
            let_wrong_pos += 1
    
    return (num_correct_pos, num_wrong_pos, num_too_small, let_correct_pos, let_wrong_pos)

def validate_solution(solution):
    # Check against all given conditions
    guesses = [
        ("35WJ", (1,0,1,0,0)),  # one num correct pos, one num too small, both letters wrong
        ("67VX", (0,0,0,0,0)),  # both nums wrong, both letters wrong
        ("37QK", (0,0,0,1,0)),  # both nums wrong, one letter correct pos
        ("75GT", (1,0,0,0,0)),  # one num correct pos, one num wrong
        ("34QE", (0,0,2,0,0)),  # both nums too small
        ("27CU", (0,0,0,0,0)),  # both nums wrong, both letters wrong
        ("56FA", (0,1,0,0,0)),  # one num wrong pos, letters too early
        ("84YT", (1,0,1,1,0))   # one num correct pos, one num too small, one letter correct pos
    ]
    
    for guess, expected in guesses:
        if check_guess(guess, solution) != expected:
            return False
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Try all possible combinations
for nums in permutations(numbers, 2):
    for lets in permutations(letters, 2):
        solution = nums[0] + nums[1] + lets[0] + lets[1]
        if validate_solution(solution):
            print(f"Found solution: {solution}")