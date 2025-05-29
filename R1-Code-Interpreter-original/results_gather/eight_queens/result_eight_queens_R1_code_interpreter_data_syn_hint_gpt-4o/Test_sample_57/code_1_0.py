def check_diagonal_constraint(positions):
    for i in range(len(positions)):
        r1, c1 = map(int, positions[i].split())
        for j in range(i + 1, len(positions)):
            r2, c2 = map(int, positions[j].split())
            if abs(r1 - r2) == abs(c1 - c2):
                return False
    return True

positions = ['0 0', '1 4', '2 7', '3 5', '4 2', '5 6', '6 1', '7 3']
is_valid = check_diagonal_constraint(positions)
print(is_valid)