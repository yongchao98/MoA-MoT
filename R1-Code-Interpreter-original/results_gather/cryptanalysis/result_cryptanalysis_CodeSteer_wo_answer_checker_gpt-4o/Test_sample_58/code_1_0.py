from itertools import permutations

# Define the initial domains
numbers = set(range(10))
letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Apply constraints from each guess
# Guess 1: 10PR
numbers -= {1, 0}
letters -= {'P', 'R'}

# Guess 2: 79EN
numbers -= {9}
letters -= {'E', 'N'}

# Guess 3: 70GE
numbers -= {0}
letters -= {'G', 'E'}

# Guess 4: 94QE
numbers -= {9, 4}
letters -= {'Q', 'E'}

# Guess 5: 05LE
numbers -= {0}
letters -= {'L', 'E'}

# Guess 6: 83CS
numbers -= {8, 3}
letters -= {'C', 'S'}

# Guess 7: 03ZF
numbers -= {0, 3}
letters -= {'F'}

# From guess 7, Z is correct and in the correct position
correct_letter = 'Z'

# Now, we need to find the correct numbers and the other letter
# We know the numbers are 6 and 7 from the analysis
possible_numbers = [6, 7]

# The remaining letters are those not eliminated and after 'Z'
possible_letters = [l for l in letters if l > correct_letter]

# Search for the correct combination
for num_perm in permutations(possible_numbers, 2):
    for letter in possible_letters:
        if letter != correct_letter:
            # Construct the password
            password = [str(num_perm[0]), str(num_perm[1]), correct_letter, letter]
            print(f"<<< {password} >>>")
            break