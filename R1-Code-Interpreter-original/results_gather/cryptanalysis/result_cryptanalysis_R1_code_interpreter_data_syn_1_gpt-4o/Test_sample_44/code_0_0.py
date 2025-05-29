# Possible numbers and letters based on the feedback
possible_numbers = [7, 6]
possible_letters = ['K', 'R', 'G']

# We know from guess 6 that one letter is correct but in the wrong position
# From guess 12, we know one letter is correct and in the correct position

# Let's deduce the correct letters
for letter1 in possible_letters:
    for letter2 in possible_letters:
        if letter1 != letter2:
            # Check against the feedback
            if (letter1 == 'K' and letter2 == 'R') or (letter1 == 'R' and letter2 == 'K'):
                # This combination satisfies the feedback from guess 6
                print([7, 6, letter1, letter2])