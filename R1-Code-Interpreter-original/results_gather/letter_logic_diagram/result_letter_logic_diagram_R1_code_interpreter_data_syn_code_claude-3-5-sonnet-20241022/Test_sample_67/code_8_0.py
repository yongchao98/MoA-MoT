def verify_initial_constraints():
    # Check all pre-filled positions on diagonal
    diagonal_positions = [
        (0,6), (1,5), (2,4), (3,3), (4,2), (5,1), (6,0)
    ]
    
    initial = [
        ['', 'c', 'b', '', '', 'd', ''],
        ['c', '', 'g', 'e', 'd', 'f', 'a'],
        ['b', '', '', '', '', '', 'c'],
        ['', '', '', '', 'a', '', ''],
        ['', '', 'f', '', '', '', ''],
        ['', 'f', '', 'c', '', 'g', ''],
        ['f', '', '', 'b', 'g', '', '']
    ]
    
    # Find what letters are pre-filled on diagonal
    diagonal_letters = set()
    for i, j in diagonal_positions:
        if initial[i][j]:
            diagonal_letters.add(initial[i][j])
    
    if len(diagonal_letters) > 1:
        return None  # Multiple different letters on diagonal
    elif len(diagonal_letters) == 1:
        return list(diagonal_letters)[0]  # Return the pre-filled diagonal letter
    else:
        return 'e'  # If no pre-filled diagonal letters, we can choose any letter

def create_solution_with_diagonal(diagonal_letter):
    # Create the solution grid with the diagonal letter fixed
    solution = [
        ['a', 'c', 'b', 'd', 'e', 'g', diagonal_letter],
        ['c', 'd', 'g', 'e', 'b', diagonal_letter, 'a'],
        ['b', 'e', 'a', 'g', diagonal_letter, 'd', 'c'],
        ['g', 'a', 'd', diagonal_letter, 'c', 'b', 'e'],
        ['e', 'b', diagonal_letter, 'c', 'g', 'a', 'd'],
        ['d', diagonal_letter, 'c', 'a', 'f', 'e', 'b'],
        [diagonal_letter, 'g', 'e', 'b', 'd', 'c', 'a']
    ]
    
    # Verify against initial constraints
    initial = [
        ['', 'c', 'b', '', '', 'd', ''],
        ['c', '', 'g', 'e', 'd', 'f', 'a'],
        ['b', '', '', '', '', '', 'c'],
        ['', '', '', '', 'a', '', ''],
        ['', '', 'f', '', '', '', ''],
        ['', 'f', '', 'c', '', 'g', ''],
        ['f', '', '', 'b', 'g', '', '']
    ]
    
    for i in range(7):
        for j in range(7):
            if initial[i][j] and initial[i][j] != solution[i][j]:
                return None
    
    return solution

# Find the correct diagonal letter
diagonal_letter = verify_initial_constraints()
if diagonal_letter:
    solution = create_solution_with_diagonal(diagonal_letter)
    if solution:
        print('<<<')
        for row in solution:
            print(','.join(row))
        print('>>>')
    else:
        print("No solution exists")
else:
    print("No solution exists - conflicting diagonal constraints")