def analyze_specific_cases(pwd):
    # Guess 1: 41CK
    def check_guess1():
        # For 41CK: one letter correct but wrong position
        # If K is in wrong position, it means K in original was in position 4 (last position)
        # So K cannot be in last position in our solution
        if pwd[3] == 'K':
            return False
        return True

    # Guess 2: 64DE
    def check_guess2():
        # Both letters must be after D and E
        # This is satisfied by both K and O
        return True

    # Guess 3: 87JY
    def check_guess3():
        # Both solutions satisfy this
        return True

    # Guess 4: 12OD
    def check_guess4():
        # For 12OD: one letter correct but wrong position
        # If O is in wrong position, it means O in original was in position 3
        # So O cannot be in position 3 in our solution
        if pwd[2] == 'O':
            return False
        return True

    return all([check_guess1(), check_guess2(), check_guess3(), check_guess4()])

solutions = []
candidates = [
    ['4', '6', 'K', 'O'],
    ['4', '6', 'O', 'K']
]

for pwd in candidates:
    if analyze_specific_cases(pwd):
        solutions.append(pwd)

print(solutions)