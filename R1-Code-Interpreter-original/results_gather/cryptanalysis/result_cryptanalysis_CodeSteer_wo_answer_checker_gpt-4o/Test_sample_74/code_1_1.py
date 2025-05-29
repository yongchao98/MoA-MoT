def is_valid_combination(numbers, letters, combination):
    num1, num2, letter1, letter2 = combination
    
    # Feedback conditions
    conditions = [
        (lambda n1, n2, l1, l2: n1 not in [8, 0] and n2 not in [8, 0] and l1 not in ['Y', 'C'] and l2 not in ['Y', 'C']),
        (lambda n1, n2, l1, l2: (n1 == 5 and n2 != 8) or (n1 != 5 and n2 == 5) and l1 not in ['L', 'T'] and l2 not in ['L', 'T']),
        (lambda n1, n2, l1, l2: n1 not in [1, 3] and n2 not in [1, 3] and l1 not in ['Q', 'W'] and l2 not in ['Q', 'W']),
        (lambda n1, n2, l1, l2: n1 not in [3, 6] and n2 not in [3, 6] and l1 not in ['I', 'L'] and l2 not in ['I', 'L']),
        (lambda n1, n2, l1, l2: (n1 == 5 and n2 != 7) or (n1 != 5 and n2 == 5) and l1 == 'F' and l2 not in ['C']),
        (lambda n1, n2, l1, l2: n1 not in [4, 6] and n2 not in [4, 6] and l1 not in ['K', 'W'] and l2 not in ['K', 'W']),
        (lambda n1, n2, l1, l2: (n1 == 5 and n2 != 0) or (n1 != 5 and n2 == 5) and l1 not in ['G', 'D'] and l2 not in ['G', 'D']),
        (lambda n1, n2, l1, l2: n1 == 5 and n2 not in [0] and l1 not in ['L', 'X'] and l2 not in ['L', 'X']),
        (lambda n1, n2, l1, l2: n1 not in [9, 6] and n2 not in [9, 6] and l1 not in ['S', 'G'] and l2 not in ['S', 'G']),
        (lambda n1, n2, l1, l2: n1 == 5 and n2 not in [4] and l1 not in ['X', 'I'] and l2 not in ['X', 'I']),
        (lambda n1, n2, l1, l2: n1 not in [6, 4] and n2 not in [6, 4] and l1 == 'F' and l2 not in ['D']),
        (lambda n1, n2, l1, l2: n1 not in [3, 1] and n2 not in [3, 1] and l1 not in ['S', 'E'] and l2 not in ['S', 'E']),
        (lambda n1, n2, l1, l2: n1 not in [0, 3] and n2 not in [0, 3] and l1 == 'Z' and l2 not in ['L'])
    ]
    
    return all(cond(num1, num2, letter1, letter2) for cond in conditions)

def find_password():
    possible_numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    possible_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 == letter2:
                        continue
                    combination = (num1, num2, letter1, letter2)
                    if is_valid_combination(possible_numbers, possible_letters, combination):
                        return combination

password = find_password()
print(f"<<< {list(password)} >>>")