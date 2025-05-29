def check_guess(guess, actual):
    # Convert numbers to integers for comparison
    g_nums = [int(guess[0]), int(guess[1])]
    a_nums = [int(actual[0]), int(actual[1])]
    g_letters = [guess[2], guess[3]]
    a_letters = [actual[2], actual[3]]
    
    # Count for numbers
    correct_pos_nums = sum(1 for i in range(2) if g_nums[i] == a_nums[i])
    correct_nums = sum(1 for i in range(2) if g_nums[i] in a_nums)
    too_large_nums = sum(1 for i in range(2) if g_nums[i] > a_nums[i])
    
    # Count for letters
    correct_pos_letters = sum(1 for i in range(2) if g_letters[i] == a_letters[i])
    correct_letters = sum(1 for i in range(2) if g_letters[i] in a_letters)
    too_early_letters = sum(1 for i in range(2) if g_letters[i] < a_letters[i])
    
    return [correct_nums, correct_pos_nums, too_large_nums, 
            correct_letters, correct_pos_letters, too_early_letters]

def verify_guess(candidate, guess, expected):
    result = check_guess(guess, candidate)
    return result == expected

def find_solution():
    # All possible numbers 0-9
    numbers = range(10)
    # All possible letters A-Z
    letters = [chr(i) for i in range(65, 91)]
    
    # Known guesses and their feedback
    conditions = [
        ("15YA", [1,0,1, 1,1,1]),
        ("07JZ", [0,0,0, 0,0,0]),
        ("70RC", [0,0,0, 0,0,2]),
        ("70WL", [0,0,0, 0,0,0]),
        ("54XP", [0,0,2, 0,0,0]),
        ("59IT", [0,0,2, 1,1,1]),
        ("41WR", [1,1,1, 0,0,0]),
        ("60RF", [0,0,0, 0,0,2]),
        ("63OA", [0,0,2, 0,0,2])
    ]
    
    # Try all possible combinations
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:  # Numbers can't repeat
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:  # Letters can't repeat
                        continue
                    
                    candidate = f"{n1}{n2}{l1}{l2}"
                    valid = True
                    
                    # Check against all conditions
                    for guess, expected in conditions:
                        if not verify_guess(candidate, guess, expected):
                            valid = False
                            break
                    
                    if valid:
                        print(f"Found solution: {candidate}")
                        return [str(n1), str(n2), l1, l2]
    
    return None

solution = find_solution()
if solution:
    print(solution)