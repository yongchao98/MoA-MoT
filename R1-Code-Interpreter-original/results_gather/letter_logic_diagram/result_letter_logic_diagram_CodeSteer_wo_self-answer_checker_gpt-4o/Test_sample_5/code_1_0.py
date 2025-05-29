def solve_puzzle(grid):
    from copy import deepcopy
    from collections import defaultdict

    # Initialize domains for each cell
    domains = defaultdict(set)
    all_letters = set('abcdefg')

    # Fill the diagonal with the same letter
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_indices if grid[i][j] != '')
    diagonal_letter = (all_letters - diagonal_letters).pop()
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Initialize domains
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                domains[(i, j)] = all_letters - row_letters - col_letters

    # Function to enforce arc consistency
    def enforce_arc_consistency():
        changes = True
        while changes:
            changes = False
            for (i, j) in list(domains.keys()):
                if len(domains[(i, j)]) == 1:
                    letter = next(iter(domains[(i, j)]))
                    for k in range(7):
                        if k != j and letter in domains[(i, k)]:
                            domains[(i, k)].remove(letter)
                            changes = True
                        if k != i and letter in domains[(k, j)]:
                            domains[(k, j)].remove(letter)
                            changes = True

    # Backtracking with constraint propagation
    def backtrack():
        if all(len(domain) == 1 for domain in domains.values()):
            return True
        # Select the variable with the smallest domain
        i, j = min((cell for cell in domains if len(domains[cell]) > 1), key=lambda x: len(domains[x]))
        for letter in domains[(i, j)].copy():
            # Try assigning the letter
            grid[i][j] = letter
            # Save the current state of domains
            saved_domains = deepcopy(domains)
            # Update domains
            domains[(i, j)] = {letter}
            enforce_arc_consistency()
            if backtrack():
                return True
            # Restore the state if failed
            grid[i][j] = ''
            domains.update(saved_domains)
        return False

    # Start solving
    enforce_arc_consistency()
    if backtrack():
        # Print the solution
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Initial grid with empty cells as ''
grid = [
    ['', 'f', 'b', 'e', 'g', 'd', 'a'],
    ['b', 'e', 'g', '', 'a', 'c', ''],
    ['', '', '', 'd', 'a', '', ''],
    ['e', '', '', 'a', 'c', 'f', ''],
    ['g', '', 'a', 'c', 'f', 'b', ''],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', '', '', 'e', 'g', 'd']
]

solve_puzzle(grid)