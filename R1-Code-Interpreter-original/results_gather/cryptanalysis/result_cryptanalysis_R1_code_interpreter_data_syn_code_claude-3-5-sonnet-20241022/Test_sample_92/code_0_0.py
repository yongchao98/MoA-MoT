def check_guess(guess, numbers_feedback, letters_feedback):
    return {
        'guess': guess,
        'numbers': numbers_feedback,
        'letters': letters_feedback
    }

guesses = [
    check_guess('58TD', 'both incorrect', 'both incorrect'),
    check_guess('52SK', 'one correct position, one incorrect', 'both incorrect'),
    check_guess('19UI', 'both incorrect', 'both incorrect'),
    check_guess('75AC', 'both incorrect', 'both too early'),
    check_guess('86WD', 'one wrong position, one too large', 'both incorrect'),
    check_guess('72KL', 'one correct position, one too large', 'both too early'),
    check_guess('08LP', 'both incorrect', 'both incorrect'),
    check_guess('51PB', 'both incorrect', 'both incorrect'),
    check_guess('96FX', 'one wrong position, one too large', 'both incorrect'),
    check_guess('02HU', 'one correct position, one too small', 'both incorrect'),
    check_guess('96FJ', 'one wrong position, one too large', 'both too early'),
    check_guess('12KC', 'one correct position, one too small', 'both too early'),
    check_guess('92TG', 'one correct position, one too large', 'both incorrect'),
    check_guess('98OJ', 'both too large', 'both incorrect'),
    check_guess('06ML', 'one wrong position, one too small', 'both too early'),
    check_guess('02VX', 'one correct position, one too small', 'one wrong position, one too late')
]

# Let's analyze the number patterns
def analyze_numbers():
    print("Number Analysis:")
    # Looking at cases where one number is correct and in correct position
    for g in guesses:
        if 'one correct position' in g['numbers']:
            print(f"Guess {g['guess']}: {g['numbers']}")

# Let's analyze the letter patterns
def analyze_letters():
    print("\nLetter Analysis:")
    # Looking at cases with useful letter information
    for g in guesses:
        if g['letters'] != 'both incorrect':
            print(f"Guess {g['guess']}: {g['letters']}")

analyze_numbers()
analyze_letters()