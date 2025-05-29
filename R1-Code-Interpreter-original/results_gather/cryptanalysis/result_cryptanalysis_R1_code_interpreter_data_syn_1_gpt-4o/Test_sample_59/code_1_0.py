# Possible numbers and letters
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

# Conditions
conditions = [
    ("20IX", (0, 0)),
    ("09YT", (1, 0)),
    ("78LK", (0, 0)),
    ("15LD", (1, 0)),
    ("03CO", (0, 1)),
    ("21XS", (1, 0)),
    ("48TP", (0, 0)),
    ("25LH", (0, 0)),
    ("72NO", (0, 1)),
    ("65RX", (0, 0))
]

# Function to check if a guess satisfies the conditions
def check_guess(guess, conditions):
    for cond, feedback in conditions:
        num_correct = sum(1 for a, b in zip(guess, cond) if a == b)
        num_incorrect = sum(1 for a in guess if a in cond) - num_correct
        if (num_correct, num_incorrect) != feedback:
            return False
    return True

# Find the correct combination
for num1 in numbers:
    for num2 in numbers:
        if num1 == num2:
            continue
        for letter1 in letters:
            for letter2 in letters:
                if letter1 == letter2:
                    continue
                guess = num1 + num2 + letter1 + letter2
                if check_guess(guess, conditions):
                    print([num1, num2, letter1, letter2])