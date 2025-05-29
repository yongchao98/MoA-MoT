def is_valid_combination(numbers, letters, combination):
    # Check each condition
    # Condition 1: 98SO
    if not (combination[0] < 9 and combination[1] < 8 and
            ((combination[2] == 'S' and combination[3] != 'O') or
             (combination[2] != 'S' and combination[3] == 'O'))):
        return False

    # Condition 2: 34ZW
    if not ((combination[0] == 3 or combination[1] == 3) and
            (combination[0] != 4 and combination[1] != 4) and
            combination[2] < 'Z' and combination[3] < 'W'):
        return False

    # Condition 3: 91HL
    if not ((combination[2] == 'H' and combination[3] != 'L') or
            (combination[2] != 'H' and combination[3] == 'L')):
        return False

    # Condition 4: 47OB
    if not ((combination[0] == 4 and combination[1] != 7) or
            (combination[0] != 4 and combination[1] == 7)):
        return False

    # Condition 5: 65NA
    if not ((combination[0] == 6 and combination[1] != 5) or
            (combination[0] != 6 and combination[1] == 5)):
        return False

    return True

def find_password():
    numbers = [2, 3, 4, 5, 6, 7]
    letters = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'S']

    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    combination = [n1, n2, l1, l2]
                    if is_valid_combination(numbers, letters, combination):
                        return combination

password = find_password()
print(f"<<< {password} >>>")