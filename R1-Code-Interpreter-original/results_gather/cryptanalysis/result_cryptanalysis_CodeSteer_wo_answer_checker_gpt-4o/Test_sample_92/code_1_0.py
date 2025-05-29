from itertools import permutations

# Define the constraints based on the feedback
constraints = [
    ("58TD", (0, 0, 0, 0)),
    ("52SK", (1, 0, 0, 0)),
    ("19UI", (0, 0, 0, 0)),
    ("75AC", (0, 0, 0, 0)),
    ("86WD", (0, 1, 0, 0)),
    ("72KL", (1, 0, 0, 0)),
    ("08LP", (0, 0, 0, 0)),
    ("51PB", (0, 0, 0, 0)),
    ("96FX", (0, 1, 0, 0)),
    ("02HU", (1, 0, 0, 0)),
    ("96FJ", (0, 1, 0, 0)),
    ("12KC", (1, 0, 0, 0)),
    ("92TG", (1, 0, 0, 0)),
    ("98OJ", (0, 0, 0, 0)),
    ("06ML", (0, 1, 0, 0)),
    ("02VX", (1, 0, 1, 0))
]

# Possible numbers and letters
numbers = "0123456789"
letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Function to check if a guess satisfies the constraints
def satisfies_constraints(guess, feedback):
    for (test_guess, (num_correct, num_wrong_pos, let_correct, let_wrong_pos)) in constraints:
        num_correct_count = sum(1 for i in range(2) if guess[i] == test_guess[i])
        num_wrong_pos_count = sum(1 for i in range(2) if guess[i] in test_guess and guess[i] != test_guess[i])
        let_correct_count = sum(1 for i in range(2, 4) if guess[i] == test_guess[i])
        let_wrong_pos_count = sum(1 for i in range(2, 4) if guess[i] in test_guess and guess[i] != test_guess[i])
        
        if (num_correct_count != num_correct or num_wrong_pos_count != num_wrong_pos or
            let_correct_count != let_correct or let_wrong_pos_count != let_wrong_pos):
            return False
    return True

# Try all permutations of two numbers and two letters
for num_perm in permutations(numbers, 2):
    for let_perm in permutations(letters, 2):
        guess = num_perm + let_perm
        if satisfies_constraints(guess, constraints):
            print(f"<<< {list(guess)} >>>")
            break