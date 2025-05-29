def flip_pieces(grid, row, col, player):
    # Directions: horizontal, vertical, and diagonal
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    opponent = 1 if player == 0 else 0

    for dr, dc in directions:
        r, c = row + dr, col + dc
        pieces_to_flip = []

        while 0 <= r < 4 and 0 <= c < 4 and grid[r][c] == opponent:
            pieces_to_flip.append((r, c))
            r += dr
            c += dc

        if 0 <= r < 4 and 0 <= c < 4 and grid[r][c] == player:
            for rr, cc in pieces_to_flip:
                grid[rr][cc] = player

# Initial grid setup
grid = [
    [1, 0, '*', '*'],
    ['*', 1, 0, '*'],
    ['*', 0, 1, '*'],
    ['*', '*', '*', '*']
]

# Round 1: Black's move at (1,2)
grid[0][1] = 0
flip_pieces(grid, 0, 1, 0)

# Round 2: White's move at (1,1)
grid[0][0] = 1
flip_pieces(grid, 0, 0, 1)

# Convert the grid to a single string
result = ','.join(str(cell) for row in grid for cell in row)
print(result)