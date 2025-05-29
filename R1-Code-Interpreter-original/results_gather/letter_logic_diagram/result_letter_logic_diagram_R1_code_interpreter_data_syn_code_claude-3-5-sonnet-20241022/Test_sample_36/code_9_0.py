def rotate_row_to_g_diagonal(row, g_position):
    """Rotate a row until 'g' is in the specified position"""
    while row[g_position] != 'g':
        row = row[1:] + [row[0]]
    return row

def create_solution():
    # Start with empty grid
    grid = [[''] * 7 for _ in range(7)]
    
    # Fill in the known pre-filled positions
    prefilled = [
        ['', '', '', 'e', 'd', 'c', ''],
        ['f', '', 'e', '', '', '', ''],
        ['', '', 'd', '', 'g', '', ''],
        ['e', '', 'c', '', 'b', 'f', ''],
        ['d', '', '', 'b', '', 'a', 'e'],
        ['c', 'g', 'b', '', 'a', 'e', 'd'],
        ['', 'b', 'f', '', '', 'd', '']
    ]
    
    # First, place 'g' on the minor diagonal
    for i in range(7):
        grid[i][6-i] = 'g'
    
    # Create valid rows that include 'g' in the correct position
    rows = [
        ['b', 'f', 'a', 'e', 'd', 'c', 'g'],  # g must be at position 6
        ['f', 'c', 'e', 'd', 'b', 'g', 'a'],  # g must be at position 5
        ['a', 'd', 'b', 'c', 'g', 'e', 'f'],  # g must be at position 4
        ['e', 'b', 'c', 'g', 'f', 'a', 'd'],  # g must be at position 3
        ['d', 'c', 'g', 'f', 'b', 'e', 'a'],  # g must be at position 2
        ['c', 'g', 'b', 'd', 'a', 'f', 'e'],  # g must be at position 1
        ['g', 'b', 'f', 'e', 'c', 'd', 'a']   # g must be at position 0
    ]
    
    # Fill the grid with these rows
    for i in range(7):
        grid[i] = rows[i]
    
    return grid

# Create and print solution
solution = create_solution()
print("<<<")
result = ""
for row in solution:
    result += ','.join(row) + '\n'
print(result.strip())
print(">>>")