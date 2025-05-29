from itertools import permutations

def is_valid_combination(combination):
    num1, num2, letter1, letter2 = combination
    
    # Check each condition
    conditions = [
        # Guess 1: 79EQ
        lambda n1, n2, l1, l2: n1 not in "79" and n2 not in "79" and l1 not in "EQ" and l2 not in "EQ",
        # Guess 2: 37TK
        lambda n1, n2, l1, l2: n1 not in "37" and n2 not in "37" and l1 < "T" and l2 < "K",
        # Guess 3: 25EO
        lambda n1, n2, l1, l2: n1 not in "25" and n2 not in "25" and l1 not in "EO" and l2 not in "EO",
        # Guess 4: 86LZ
        lambda n1, n2, l1, l2: (n1 == "6" or n2 == "6") and n1 not in "86" and n2 not in "86" and l1 < "L" and l2 < "Z",
        # Guess 5: 43WQ
        lambda n1, n2, l1, l2: n1 not in "43" and n2 not in "43" and l1 < "W" and l2 < "Q",
        # Guess 6: 31AO
        lambda n1, n2, l1, l2: n1 not in "31" and n2 not in "31" and (l1 == "A" or l2 == "A") and l1 < "O" and l2 < "O",
        # Guess 7: 89GW
        lambda n1, n2, l1, l2: n1 not in "89" and n2 not in "89" and (l1 == "G" or l2 == "G") and l1 < "W" and l2 < "W"
    ]
    
    return all(cond(num1, num2, letter1, letter2) for cond in conditions)

def find_password():
    possible_numbers = "0123456789"
    possible_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # Filter based on initial conditions
    possible_numbers = [n for n in possible_numbers if n not in "379258"]
    possible_letters = [l for l in possible_letters if l < "T" and l not in "EQOZKW"]
    
    # Generate all permutations of two numbers and two letters
    for num_comb in permutations(possible_numbers, 2):
        for letter_comb in permutations(possible_letters, 2):
            combination = num_comb + letter_comb
            if is_valid_combination(combination):
                return combination

    return "No valid combination found"

# Find and print the password
password = find_password()
print(f"<<< {list(password)} >>>")