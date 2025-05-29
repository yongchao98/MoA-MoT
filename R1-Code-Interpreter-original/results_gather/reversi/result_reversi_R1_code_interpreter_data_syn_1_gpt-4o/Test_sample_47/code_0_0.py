def flip_pieces(grid, row, col, player):
    n = len(grid)
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    
    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []
        
        while 0 <= r < n and 0 <= c < n and grid[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc
        
        if 0 <= r < n and 0 <= c < n and grid[r][c] == player:
            for rr, cc in pieces_to_flip:
                grid[rr][cc] = player

def apply_move(grid, row, col, player):
    grid[row][col] = player
    flip_pieces(grid, row, col, player)

# Initial grid configuration
grid = [
    ['*', '*', '1', '*'],
    ['0', '0', '1', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '*', '0']
]

# Apply moves
apply_move(grid, 1, 0, '0')  # Round 1: Black at (2,1)
apply_move(grid, 0, 2, '1')  # Round 2: White at (1,3)
apply_move(grid, 3, 3, '0')  # Round 3: Black at (4,4)

# Flatten the grid for the final output
final_grid = [cell for row in grid for cell in row]
print(final_grid)