from itertools import product

def check_guess(guess, test):
    # Returns True if the test combination satisfies all conditions for the given guess and feedback
    
    # Helper function to check number conditions
    def check_numbers(guess_nums, test_nums, feedback):
        if feedback == "both incorrect":
            return test_nums[0] != guess_nums[0] and test_nums[1] != guess_nums[1]
        elif feedback == "both incorrect and too large":
            return all(t != g and t < g for t, g in zip(test_nums, guess_nums))
        elif feedback == "one correct in position":
            return (test_nums[0] == guess_nums[0] and test_nums[1] != guess_nums[1]) or \
                   (test_nums[1] == guess_nums[1] and test_nums[0] != guess_nums[0])
        elif feedback == "one correct wrong position":
            return (test_nums[0] == guess_nums[1] or test_nums[1] == guess_nums[0]) and \
                   test_nums[0] != guess_nums[0] and test_nums[1] != guess_nums[1]
    
    # Helper function to check letter conditions
    def check_letters(guess_lets, test_lets, feedback):
        if feedback == "both incorrect":
            return test_lets[0] != guess_lets[0] and test_lets[1] != guess_lets[1]
        elif feedback == "one correct in position":
            return (test_lets[0] == guess_lets[0] and test_lets[1] != guess_lets[1]) or \
                   (test_lets[1] == guess_lets[1] and test_lets[0] != guess_lets[0])
        elif feedback == "one correct wrong position":
            return (test_lets[0] == guess_lets[1] or test_lets[1] == guess_lets[0]) and \
                   test_lets[0] != guess_lets[0] and test_lets[1] != guess_lets[1]

    # Check all conditions
    conditions = [
        ("86AF", ("both incorrect", "both incorrect")),
        ("98LF", ("both incorrect and too large", "both incorrect")),
        ("20XK", ("one correct in position", "both incorrect")),
        ("38TM", ("both incorrect", "both incorrect")),
        ("81JK", ("both incorrect", "both incorrect")),
        ("23VQ", ("one correct in position", "both incorrect")),
        ("85IC", ("both incorrect", "both incorrect")),
        ("03OE", ("both incorrect", "both incorrect")),
        ("78WK", ("one correct wrong position", "one correct in position")),
        ("68KC", ("both incorrect", "both incorrect")),
        ("47OM", ("one correct in position", "both incorrect")),
        ("58ZQ", ("both incorrect", "both incorrect")),
        ("49DW", ("both incorrect", "one correct wrong position"))
    ]

    for guess_str, (num_feedback, let_feedback) in conditions:
        guess_nums = [int(guess_str[0]), int(guess_str[1])]
        guess_lets = [guess_str[2], guess_str[3]]
        
        if not check_numbers(guess_nums, [test[0], test[1]], num_feedback):
            return False
        if not check_letters(guess_lets, [test[2], test[3]], let_feedback):
            return False
    
    return True

# Generate all possible combinations
numbers = range(10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

valid_solutions = []
for n1, n2 in product(numbers, numbers):
    if n1 == n2:  # Skip if numbers repeat
        continue
    for l1, l2 in product(letters, letters):
        if l1 == l2:  # Skip if letters repeat
            continue
        test = [n1, n2, l1, l2]
        if check_guess(None, test):
            valid_solutions.append(test)

print(f"Found {len(valid_solutions)} solution(s):")
for solution in valid_solutions:
    print(solution)