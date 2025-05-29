def flip_pieces(grid, row, col, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
    n = len(grid)
    
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

def play_game():
    grid = [
        ['*', '*', '*', '*'],
        ['*', '1', '1', '1'],
        ['*', '0', '0', '0'],
        ['*', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black, 0) at (3, 4)
    grid[2][3] = '0'
    flip_pieces(grid, 2, 3, '0')
    
    # Round 2: Opponent's move (White, 1) at (2, 4)
    grid[1][3] = '1'
    flip_pieces(grid, 1, 3, '1')
    
    # Flatten the grid for the final output
    flat_grid = [cell for row in grid for cell in row]
    return flat_grid

result = play_game()
print(result)