def check_guess(guess, actual):
    # Returns feedback similar to the given conditions
    correct_nums = 0
    correct_pos_nums = 0
    too_large_nums = 0
    correct_letters = 0
    correct_pos_letters = 0
    too_early_letters = 0
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            correct_nums += 1
            correct_pos_nums += 1
        elif guess[i] in actual[:2]:
            correct_nums += 1
        elif int(guess[i]) > int(actual[i]):
            too_large_nums += 1
            
    # Check letters
    for i in range(2,4):
        if guess[i] == actual[i]:
            correct_letters += 1
            correct_pos_letters += 1
        elif guess[i] in actual[2:]:
            correct_letters += 1
        elif guess[i] < actual[i]:
            too_early_letters += 1
            
    return (correct_nums, correct_pos_nums, too_large_nums, 
            correct_letters, correct_pos_letters, too_early_letters)

def matches_conditions(candidate):
    # Test against all known guesses
    guesses = [
        ("15YA", (1,0,1, 1,1,1)),
        ("07JZ", (0,0,0, 0,0,0)),
        ("70RC", (0,0,0, 0,0,2)),
        ("70WL", (0,0,0, 0,0,0)),
        ("54XP", (0,0,2, 0,0,0)),
        ("59IT", (0,0,2, 1,1,1)),
        ("41WR", (1,1,1, 0,0,0)),
        ("60RF", (0,0,0, 0,0,2)),
        ("63OA", (0,0,2, 0,0,2))
    ]
    
    for guess, expected in guesses:
        if check_guess(guess, candidate) != expected:
            return False
    return True

# Generate all possible combinations
possible_numbers = [str(i).zfill(2) for i in range(100)]
possible_letters = [chr(i) for i in range(65, 91)]  # A-Z
solutions = []

for nums in possible_numbers:
    if nums[0] == nums[1]:  # Skip repeating numbers
        continue
    for l1 in possible_letters:
        for l2 in possible_letters:
            if l1 == l2:  # Skip repeating letters
                continue
            candidate = nums + l1 + l2
            if matches_conditions(candidate):
                solutions.append(candidate)

print(solutions)