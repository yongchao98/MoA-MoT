def flip_pieces(grid, player, row, col):
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
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
        ['*', '0', '*', '*'],
        ['*', '0', '0', '*'],
        ['1', '0', '1', '*'],
        ['0', '*', '*', '*']
    ]
    
    # Round 1: Your move (Black) at (1, 2)
    grid[0][1] = '0'
    flip_pieces(grid, '0', 0, 1)
    
    # Round 2: Opponent's move (White) at (3, 1)
    grid[2][0] = '1'
    flip_pieces(grid, '1', 2, 0)
    
    # Round 3: Your move (Black) at (4, 1)
    grid[3][0] = '0'
    flip_pieces(grid, '0', 3, 0)
    
    # Flatten the grid to the required format
    result = [cell for row in grid for cell in row]
    print(result)

play_game()