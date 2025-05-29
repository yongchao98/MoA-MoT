def solve_puzzle():
    # Initial grid with given letters
    grid = [
        ['d', '', '', 'f', 'a', '', ''],
        ['b', '', 'f', '', '', 'e', 'd'],
        ['', 'f', 'a', 'c', 'e', 'd', 'b'],
        ['f', 'a', 'c', '', '', '', ''],
        ['a', '', 'e', '', '', 'g', 'f'],
        ['c', 'e', '', 'b', '', '', 'a'],
        ['e', 'd', '', 'g', 'f', 'a', '']
    ]

    # Determine the letter for the minor diagonal
    # Check which letter is most suitable for the diagonal
    # We can try 'c' as it is not present in any of the diagonal positions
    diagonal_letter = 'c'
    
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Function to find missing letters in a list
    def find_missing_letters(lst):
        all_letters = set('abcdefg')
        return list(all_letters - set(lst))

    # Fill the grid row by row
    for i in range(7):
        missing_letters = find_missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing_letters.pop()

    # Transpose the grid to check columns
    transposed_grid = list(map(list, zip(*grid)))

    # Fill the grid column by column
    for j in range(7):
        missing_letters = find_missing_letters(transposed_grid[j])
        for i in range(7):
            if transposed_grid[j][i] == '':
                transposed_grid[j][i] = missing_letters.pop()

    # Transpose back to original grid
    grid = list(map(list, zip(*transposed_grid)))

    # Print the final grid
    for row in grid:
        print(','.join(row))

solve_puzzle()