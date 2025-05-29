from itertools import permutations

# Possible numbers and letters based on deductions
possible_numbers = [0, 4, 5, 6]
possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']

# Generate all permutations of two numbers and two letters
number_permutations = permutations(possible_numbers, 2)
letter_permutations = permutations(possible_letters, 2)

# Function to check if a combination satisfies all conditions
def check_combination(numbers, letters):
    # Convert numbers and letters to strings for comparison
    num_str = ''.join(map(str, numbers))
    let_str = ''.join(letters)
    
    # Check each condition
    # Condition 1: 79EQ
    if any(n in num_str for n in '79') or any(l in let_str for l in 'EQ'):
        return False
    
    # Condition 2: 37TK
    if any(n in num_str for n in '37') or any(l in let_str for l in 'TK'):
        return False
    
    # Condition 3: 25EO
    if any(n in num_str for n in '25') or any(l in let_str for l in 'EO'):
        return False
    
    # Condition 4: 86LZ
    if ('8' in num_str and '6' not in num_str) or ('6' in num_str and '8' not in num_str) or any(l in let_str for l in 'LZ'):
        return False
    
    # Condition 5: 43WQ
    if any(n in num_str for n in '43') or any(l in let_str for l in 'WQ'):
        return False
    
    # Condition 6: 31AO
    if any(n in num_str for n in '31') or ('A' in let_str and 'O' not in let_str) or ('O' in let_str and 'A' not in let_str):
        return False
    
    # Condition 7: 89GW
    if any(n in num_str for n in '89') or ('G' in let_str and 'W' not in let_str) or ('W' in let_str and 'G' not in let_str):
        return False
    
    return True

# Iterate through all permutations and find the correct combination
for numbers in number_permutations:
    for letters in letter_permutations:
        if check_combination(numbers, letters):
            print([str(numbers[0]), str(numbers[1]), letters[0], letters[1]])
            break