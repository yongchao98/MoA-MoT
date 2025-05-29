def print_grid(grid):
    for row in grid:
        print(''.join(row))

def flip_pieces(grid, player, opponent, row, col):
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
    n = 4
    grid = [
        ['*', '*', '1', '*'],
        ['0', '0', '1', '*'],
        ['*', '0', '0', '*'],
        ['*', '*', '*', '0']
    ]
    
    # Round 1: Your move (Black) at (2, 1)
    grid[1][0] = '0'
    flip_pieces(grid, '0', '1', 1, 0)
    
    # Round 2: Opponent's move (White) at (1, 3)
    grid[0][2] = '1'
    flip_pieces(grid, '1', '0', 0, 2)
    
    # Round 3: Your move (Black) at (4, 4)
    grid[3][3] = '0'
    flip_pieces(grid, '0', '1', 3, 3)
    
    # Print the final grid
    print_grid(grid)

play_game()