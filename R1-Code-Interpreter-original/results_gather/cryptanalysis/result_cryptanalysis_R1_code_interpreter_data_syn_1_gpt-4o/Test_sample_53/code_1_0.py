import itertools

# Possible numbers and letters
numbers = [8, 9]
letters = [chr(i) for i in range(ord('N'), ord('Z') + 1)]

# Check all permutations of numbers and letters
for num_perm in itertools.permutations(numbers, 2):
    for letter_perm in itertools.permutations(letters, 2):
        # Construct the guess
        guess = list(num_perm) + list(letter_perm)
        
        # Check against all conditions
        if (guess[0] > 3 and guess[1] > 0 and
            (guess[2] == 'A' or guess[3] == 'A') and
            (guess[2] != 'U' and guess[3] != 'U') and
            (guess[2] != 'K' and guess[3] != 'K') and
            (guess[2] != 'I' and guess[3] != 'I') and
            (guess[2] != 'Q' and guess[3] != 'Q') and
            (guess[2] != 'B' and guess[3] != 'B') and
            (guess[2] != 'H' and guess[3] != 'H') and
            (guess[2] != 'C' and guess[3] != 'C') and
            (guess[2] != 'M' and guess[3] != 'M') and
            (guess[2] != 'I' and guess[3] != 'I')):
            print(guess)