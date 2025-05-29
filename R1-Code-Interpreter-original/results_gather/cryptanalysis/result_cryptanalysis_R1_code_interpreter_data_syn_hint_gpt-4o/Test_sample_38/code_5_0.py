from itertools import permutations

# Possible numbers and letters based on the conditions
possible_numbers = [0, 2, 3]
possible_letters = ['K', 'L', 'G']

# Check each permutation of numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        # Construct the guess
        guess = list(num_perm) + list(letter_perm)
        
        # Check against the conditions
        # Condition 1: Numbers < 4, Letters < Y, Q
        if guess[0] < 4 and guess[1] < 4:
            # Condition 2: One letter is L or Q in the correct position
            if (guess[2] == 'L' or guess[3] == 'L'):
                # Condition 3: One number is 2 or 3 but in the wrong position
                if (guess[0] == 2 or guess[1] == 2 or guess[0] == 3 or guess[1] == 3):
                    # Condition 4: One number is 0 or 9 but in the wrong position
                    if (guess[0] == 0 or guess[1] == 0):
                        # Condition 4: One letter is K or G in the correct position
                        if (guess[2] == 'K' or guess[3] == 'K'):
                            print(guess)
                            break