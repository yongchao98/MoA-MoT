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
        ['*', '*', '*', '*', '*', '*'],
        ['*', '*', '0', '*', '*', '*'],
        ['*', '*', '0', '0', '*', '*'],
        ['*', '1', '0', '1', '*', '*'],
        ['*', '*', '0', '*', '*', '*'],
        ['*', '*', '*', '*', '*', '*']
    ]
    
    # Round 1: Black plays at (2,3)
    grid[1][2] = '0'
    flip_pieces(grid, 1, 2, '0')
    
    # Round 2: White plays at (4,2)
    grid[3][1] = '1'
    flip_pieces(grid, 3, 1, '1')
    
    # Round 3: Black plays at (5,3)
    grid[4][2] = '0'
    flip_pieces(grid, 4, 2, '0')
    
    # Flatten the grid to a single list
    result = [cell for row in grid for cell in row]
    print(result)

play_game()