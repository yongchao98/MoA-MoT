def is_valid_combination(num1, num2, letter1, letter2):
    # Feedback conditions
    conditions = [
        # Guess 1: 80YC - Both numbers and letters are incorrect
        lambda n1, n2, l1, l2: n1 not in [8, 0] and n2 not in [8, 0] and l1 not in ['Y', 'C'] and l2 not in ['Y', 'C'],
        # Guess 2: 58LT - One number correct but wrong position, one number too large, both letters incorrect
        lambda n1, n2, l1, l2: ((n1 == 5 and n2 != 8) or (n1 != 5 and n2 == 5)) and n1 != 8 and n2 != 8 and l1 not in ['L', 'T'] and l2 not in ['L', 'T'],
        # Guess 3: 13QW - Both numbers and letters are incorrect
        lambda n1, n2, l1, l2: n1 not in [1, 3] and n2 not in [1, 3] and l1 not in ['Q', 'W'] and l2 not in ['Q', 'W'],
        # Guess 4: 36IL - Both numbers and letters are incorrect
        lambda n1, n2, l1, l2: n1 not in [3, 6] and n2 not in [3, 6] and l1 not in ['I', 'L'] and l2 not in ['I', 'L'],
        # Guess 5: 57CF - One number correct but wrong position, one letter correct in position, one letter too early
        lambda n1, n2, l1, l2: ((n1 == 5 and n2 != 7) or (n1 != 5 and n2 == 5)) and n1 != 7 and n2 != 7 and l1 == 'F' and l2 > 'C',
        # Guess 6: 46KW - Both numbers and letters are incorrect
        lambda n1, n2, l1, l2: n1 not in [4, 6] and n2 not in [4, 6] and l1 not in ['K', 'W'] and l2 not in ['K', 'W'],
        # Guess 7: 50GD - One number correct but wrong position, one number too small, both letters incorrect
        lambda n1, n2, l1, l2: ((n1 == 5 and n2 != 0) or (n1 != 5 and n2 == 5)) and n1 != 0 and n2 != 0 and l1 not in ['G', 'D'] and l2 not in ['G', 'D'],
        # Guess 8: 05LX - One number correct in position, one number too small, both letters incorrect
        lambda n1, n2, l1, l2: n1 == 5 and n2 > 0 and l1 not in ['L', 'X'] and l2 not in ['L', 'X'],
        # Guess 9: 96SG - Both numbers too large, both letters incorrect
        lambda n1, n2, l1, l2: n1 < 9 and n2 < 6 and l1 not in ['S', 'G'] and l2 not in ['S', 'G'],
        # Guess 10: 45XI - One number correct in position, one number incorrect, both letters incorrect
        lambda n1, n2, l1, l2: n1 == 5 and n2 != 4 and l1 not in ['X', 'I'] and l2 not in ['X', 'I'],
        # Guess 11: 64DF - Both numbers incorrect, one letter correct in position, one letter too early
        lambda n1, n2, l1, l2: n1 not in [6, 4] and n2 not in [6, 4] and l1 == 'F' and l2 > 'D',
        # Guess 12: 31SE - Both numbers and letters are incorrect
        lambda n1, n2, l1, l2: n1 not in [3, 1] and n2 not in [3, 1] and l1 not in ['S', 'E'] and l2 not in ['S', 'E'],
        # Guess 13: 03LZ - Both numbers incorrect, one letter correct but wrong position
        lambda n1, n2, l1, l2: n1 not in [0, 3] and n2 not in [0, 3] and l1 == 'Z' and l2 not in ['L']
    ]
    
    for i, cond in enumerate(conditions):
        if not cond(num1, num2, letter1, letter2):
            print(f"Combination {num1}{num2}{letter1}{letter2} failed at condition {i+1}")
            return False
    return True

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
                    if is_valid_combination(num1, num2, letter1, letter2):
                        print(f"Checking combination: {num1}{num2}{letter1}{letter2}")
                        return [str(num1), str(num2), letter1, letter2]

password = find_password()
print(f"<<< {password} >>>")